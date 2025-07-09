import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg,einsum,conj,sqrt
from scipy.linalg import solve
from scipy.optimize import minimize, approx_fprime
import h5py
import sys
class RotheSolver:
    """
    This class implements the Rothe method for the time-dependent Schrödinger equation.
    Essentially, it does the following:
    1. The user has to provide a function f()
    """
    def initialize_matrices(self,S,H,H2,Sd,Hd,H2d,dt,name):
        """
        S:  Function of p, returns complex (n,n) overlap matrix.
        H:  Function of p AND t, returns complex (n,n) Hamiltonian matrix.
        H2: Function of p AND t, returns complex (n,n) Hamiltonian squared matrix.
        Sd: Derivative of S. Function of p, returns a 3-tensor of dimension (n,n,N_b)
            where the element (x,y,z) is the derivative of S(x,y) with respect to z'th nonlinear coefficient of the y'th Gaussian
            (we can obtain the derivative of the z'th nonlinear coefficient of the x'th Gaussian as (y,x,z)^*).
            In general, there are (n,n,n,N_b) derivatives, but we assume that all Gaussians are independent.
        Hd: Derivative of H.Function of p and t, returns a 3-tensor of dimension (n,n,N_b).
            Same structure as Sd, but for the Hamiltonian matrix.
        H2d: Derivative of H2. Function of p and t, returns a 3-tensor of dimension (n,n,N_b).
            Same structure as Sd, but for the Hamiltonian squared matrix.
        dt: Time step for the Rothe method.
        name: Name of the solver instance, used for saving results.
        """ 
        self.S = S
        self.H = H
        self.H2 = H2
        self.Sd = Sd
        self.Hd = Hd
        self.H2d = H2d
        self.dt = dt
        self.name = name
    def initialize_parameters(self,p,c):
        """
        p:  REAL array, the initial nonlinear coefficients.
            Shape: (n, N_b), where n is the number of Gaussians and N_b is the number of basis functions per Gaussian.
        c:  COMPLEX array, the initial linear coefficients
            Shape: (n)
        """
        self.p = p
        self.c = c
        self.n = self.p.shape[0]
        self.N_b = self.p.shape[1]
    def initialize_WF_object(self,WF):
        """
        WF: An object that contains the corresponding functions in a more "condensed" way.
        """
        self.WF=WF
        self.p=WF.p
        self.c=WF.c
        self.S=WF.S
        self.H= WF.H
        self.H2 = WF.H2
        self.Sd = WF.Sd
        self.Hd = WF.Hd
        self.H2d = WF.H2d
        self.dt = WF.dt
        self.name = WF.name
        self.n = self.p.shape[0]
        self.N_b = self.p.shape[1]
    def propagate(self,t,dt,starting_guess=None):
        """
        Propagates the system using the Rothe method with a time step of dt from time t to time t+dt.
        This method will update the coefficients p and c based on the Rothe method.
        """
        self.n = self.p.shape[0]
        self.N_b = self.p.shape[1]
        if starting_guess is None:
            p_new_init=self.p.copy() #Initial guess for the new nonlinear coefficients.
        else:
            p_new_init=starting_guess
        def error(p_new_flat):
            p_new= p_new_flat.reshape((self.n,self.N_b)) #Reshape the flat array to the original shape.
            return self.calculate_Rothe_error(t,dt,self.c,p_new,return_cnew=False)
        def gradient_numerical(p_new_flat):
            return approx_fprime(p_new_flat,error,1e-8) #Numerical gradient of the error function.
        def error_and_gradient(p_new_flat):
            p_new= p_new_flat.reshape((self.n,self.N_b)) #Reshape the flat array to the original shape.
            rothe_error, rothe_gradient = self.calculate_Rothe_error_and_gradient(t,dt,self.c,p_new)
            return rothe_error, rothe_gradient
        if self.Sd is None:
            print("I will be very slow, because I have to calculate the gradient numerically.")
            def gradient_function_numerical(p_new_flat):
                return gradient_numerical(p_new_flat)
            gradient_function = gradient_function_numerical
        else:
            def gradient_function_analytic(p_new_flat):
                p_new= p_new_flat.reshape((self.n,self.N_b)) #Reshape the flat array to the original shape.
                rothe_error, rothe_gradient = self.calculate_Rothe_error_and_gradient(t,dt,self.c,p_new)
                return rothe_gradient
            gradient_function = gradient_function_analytic
        p_new_init_flat=p_new_init.flatten() #Flatten the initial guess for the new nonlinear coefficients.
        initial_error=error(p_new_init_flat) #Initial error.
        initial_gradient=gradient_function(p_new_init_flat)
        #gradient_approx=gradient_numerical(p_new_init_flat) #Numerical gradient of the initial guess.
        #print("Gradient difference: ", np.linalg.norm(initial_gradient-gradient_approx))
        #if np.linalg.norm(initial_gradient-gradient_approx) > 1e-6:
        #    print("Warning: The analytical gradient does not match the numerical gradient. This might indicate a bug in the code.")
        #    print("Analytical gradient: ", initial_gradient)
        #    print("Numerical gradient: ", gradient_approx)
        #    sys.exit()
        epsilon = 1e-16
        hess_inv0 = np.diag(1.0 / (np.abs(np.asarray(initial_gradient)) + epsilon))
        # We use scipy.optimize.minimize to find the new coefficients.
        if self.Sd is None:
            result = minimize(error, p_new_init_flat, method='BFGS', 
                            jac=gradient_function, options={"hess_inv0": hess_inv0,'disp': False, "gtol":1e-14})
        else:
            result = minimize(error_and_gradient, p_new_init_flat, method='BFGS', 
                            jac=True, options={"hess_inv0": hess_inv0,'disp': False, "gtol":1e-14})
        p_new_flat = result.x
        p_new= p_new_flat.reshape((self.n,self.N_b)) #Reshape the flat array to the original shape.
        rothe_error, c_new = self.calculate_Rothe_error(t,dt,self.c,p_new,return_cnew=True)
        
        print("Rothe error going to time %.3f: %.2e (Initial: %.2e)" % (t+dt, rothe_error, initial_error))
        print(p_new)
        return p_new, c_new, rothe_error
    def propagate_Nsteps_and_store(
            self, t0, N, *, dt=None,
            filename=None, group_prefix="step",starting_guess_method="linear"):
        """
        Evolve the system N steps starting at t0 and stream the coefficients
        to an HDF5 file.

        Parameters
        ----------
        t0 : float
            Initial (physical) time.
        N  : int
            Number of propagation steps **after** the initial state
            (total frames written = N+1).
        dt : float, optional
            Time-step to use (default: self.dt).
        filename : str, optional
            `.h5` file to create/append.  Default = f"{self.name}.h5".
        group_prefix : str, optional
            Each frame is written into  f"/{group_prefix}_0000", etc.
        starting_guess_method : str, optional
            Method to use for the starting guess of the new nonlinear coefficients.
            Options are "linear" (use the new coefficients as starting guess for the next step)
            or "difference" (use the difference between the new and old coefficients as starting guess)
        """
        dt = self.dt if dt is None else float(dt)
        fname = f"{self.name}.h5" if filename is None else filename

        # open HDF5 file ( ‘r+’ if it already exists so we can append )
        with h5py.File(fname, "w") as f:
            # write/overwrite some global metadata once
            f.attrs.update(dict(
                solver="Rothe", dt=dt, Nb=self.N_b, created_by=self.name
            ))

            # ----------------------- write initial frame -------------
            t = t0
            g = f.create_group(f"{group_prefix}_{0:06d}")
            g.attrs["time"] = t
            g.create_dataset("p", data=self.p)                  # real
            g.create_dataset("c", data=self.c.astype(np.complex128))

            # ----------------------- propagate & stream -------------
            print("Initial state")
            print(self.p)
            starting_guess = self.p.copy()  # Initial guess for the new nonlinear coefficients.
            for step in range(1, N + 1):
                p_new, c_new, err = self.propagate(t, dt,starting_guess=starting_guess)
                if starting_guess_method == "linear":
                    starting_guess= p_new.copy()  # Use the new coefficients as starting guess for the next step.
                elif starting_guess_method == "difference":
                    starting_guess=p_new+(p_new- self.p) #Use the difference as starting guess for the next step.
                # advance internal state
                self.p, self.c = p_new, c_new
                t += dt

                # store
                g = f.create_group(f"{group_prefix}_{step:06d}")
                g.attrs["time"] = t
                g.create_dataset("p", data=p_new)
                g.create_dataset("c", data=c_new.astype(np.complex128))
                g.attrs["rothe_error"] = np.asarray(err)

            print(f"✓ wrote {N+1} frames to {fname}")

    # ────────────────────────────────────────────────────────────────
    @staticmethod
    def read_data(time, filename:str, *, atol=1e-8, group_prefix="step"):
        """
        Fetch the nonlinear (p) and linear (c) coefficients that were saved
        closest to `time` (|t_saved − time| < atol).

        Returns
        -------
        p : np.ndarray   (n, N_b)
        c : np.ndarray   (n,)   complex
        saved_time : float      the exact timestamp found in the file
        """
        with h5py.File(filename, "r") as f:
            for key in f:
                if not key.startswith(group_prefix):
                    continue
                g = f[key]
                t_saved = g.attrs["time"]
                if abs(t_saved - time) <= atol:
                    p = g["p"][...]
                    c = g["c"][...]
                    return p, c, t_saved
        raise ValueError(
            f"No frame with |t−{time}| < {atol} found in {filename}"
        )
    def calculate_Rothe_error(self,t,dt,c_old,p_new,return_cnew=False,lambda_=1e-16):
        """
        Calculate Rothe error going from time $t$ to $t+ Delta t$. 
        This function assumes that self.p and self.c represent the wave function at the PREVIOUS time step.
        """
        ngo=num_gauss_old= self.p.shape[0]
        ngn=num_gauss_new=p_new.shape[0]
        p_concat=np.concatenate((self.p,p_new),axis=0) #New shape: (2*n,N_b)
        S_full= self.S(p_concat)
        H_full = self.H(p_concat,t+dt/2) #Half time step Hamiltonian
        H2_full = self.H2(p_concat,t+dt/2)
        S_tilde_full=S_full + dt**2/4* H2_full
        rho_mat   = (S_full - (dt**2/4) * H2_full + 1j*dt*H_full)[:ngo, ngo:]
        rho_vec   = np.conj(rho_mat).T @ c_old
        S_tilde_full=0.5*(S_tilde_full + np.conj(S_tilde_full.T))
        S_reg= S_tilde_full + np.eye(S_tilde_full.shape[0])*lambda_ #Regularization term to avoid singular matrices.
        try:
            c_new = solve(S_reg[ngo:,ngo:],rho_vec,assume_a='pos')
        except np.linalg.LinAlgError as e:
            print("Error in solving the linear system: ", e)
            print(S_reg[ngo:,ngo:])
            import sys
            sys.exit(0)
        overlap_term = np.conj(c_old)@S_tilde_full[:ngo,:ngo]@c_old
        #projection_term=2*conj(rho_vec).T@c_new-conj(c_new).T@S_tilde_full[ngo:,ngo:]@c_new
        projection_term = np.real(np.conj(rho_vec).T @ c_new )
        rothe_error = overlap_term - projection_term
        rothe_error=np.real(overlap_term - projection_term)
        
        if return_cnew:
            return rothe_error, c_new
        else:
            return rothe_error
    def calculate_Rothe_error_and_gradient(self, t, dt, c_old, p_new, λ=1e-14):
        """
        Return Rothe error and its analytic gradient, using a Tikhonov
        regularisation  λ·I  on the new-new block.
        """
        ngo  = num_gauss_old = self.p.shape[0]     # old Gaussians
        Nb   = self.p.shape[1]
        ngn  = num_gauss_new = p_new.shape[0]

        # ────────────────── build full matrices ────────────────────────────────
        p_concat = np.concatenate((self.p, p_new), axis=0)
        S_full   = self.S(p_concat)
        H_full   = self.H(p_concat,  t + dt/2)
        H2_full  = self.H2(p_concat, t + dt/2)

        # enforce Hermiticity (round-off safety)
        for M in (S_full, H_full, H2_full):
            M[:] = 0.5 * (M + M.conj().T)

        Σ_full   = S_full + (dt**2/4) * H2_full
        Σ_nn     = Σ_full[ngo:, ngo:]
        Σ_oo     = Σ_full[:ngo, :ngo]

        # ────────────────── right-hand side and regularised solve ──────────────
        R_mat    = (S_full - (dt**2/4) * H2_full + 1j*dt*H_full)[:ngo, ngo:]
        ρ        = R_mat.conj().T @ c_old

        simga_tilde= Σ_nn + λ * np.eye(ngn, dtype=Σ_nn.dtype)     # ⇢ CHANGED
        c_new = solve(simga_tilde, ρ, assume_a='pos')              # ⇢ CHANGED

        # ────────────────── scalar error ---------------------------------------
        overlap  = np.real(np.vdot(c_old, Σ_oo @ c_old))
        proj     = np.real(np.vdot(ρ, c_new))
        rothe_error = overlap - proj
        # ────────────────── derivative tensors ---------------------------------
        Sd   = self.Sd(p_concat)
        Hd   = self.Hd(p_concat,  t + dt/2)
        H2d  = self.H2d(p_concat, t + dt/2)

        #   dΣ(k,a,b) = ∂(S + Δt²/4 H²) / ∂p_ab   (k row in NEW block)
        dΣ_cols = Sd[ngo:, ngo:, :] + (dt**2/4.0) * H2d[ngo:, ngo:, :]   # (ngn, ngn, Nb)
        #   dR(i,a,b) = ∂(S − Δt²/4 H² + iΔt H) / ∂p_ab,  i in OLD block
        dR_cols = ( Sd[:ngo,  ngo:, :]
                    - (dt**2/4.0) * H2d[:ngo, ngo:, :]
                    + 1j*dt * Hd[:ngo,  ngo:, :] )                      # (ngo, ngn, Nb)

        # ── helper contractions (pure matmuls, no einsum) ──────────────────────
        dot_cs = (np.conj(c_new)[None, :]                      # (1, ngn)
                @ dΣ_cols.reshape(ngn, ngn*Nb)).reshape(ngn, Nb)

        dRho   = (c_old[None, :]                               # (1, ngo)
                @ np.conj(dR_cols).reshape(ngo, ngn*Nb)
                ).reshape(ngn, Nb)

        dΣ_diag = dΣ_cols[np.arange(ngn), np.arange(ngn), :]   # (ngn, Nb)

        # ── gradient components ────────────────────────────────────────────────
        term1 = 2.0 * np.real(c_new[:, None] * dot_cs) \
                - np.real((np.abs(c_new)**2)[:, None] * dΣ_diag)

        term2 = np.real(np.conj(c_new)[:, None] * dRho)

        rothe_gradient = (term1 - 2.0*term2).ravel()
        return rothe_error, rothe_gradient
    
if __name__ == "__main__":
    from HOscillator_1D import HOscillator_1D
    from RotheSolver     import RotheSolver     # includes the I/O helpers

    # ---------- initial wave-function parameters ------------------------
    dt     = 0.01
    omega  = 1.0
    p_init = np.array([[1/np.sqrt(2), 0,1,0],[1/np.sqrt(2), 0,-1,0]
                    
                    ])          # p
    c_init = np.array([1,1])  # c

    # ---------- build WF + solver ---------------------------------------
    wf     = HOscillator_1D(p_init, c_init,
                            omega=omega, dt=dt, name="gauss_x2")
    solver = RotheSolver()
    solver.initialize_WF_object(wf)
    initial_error,initial_gradient = solver.calculate_Rothe_error_and_gradient(0,dt,c_init,p_init)
    print(solver.S(p_init))
    print(solver.H(p_init,0))
    print(solver.H2(p_init,0))
    print("Initial error: ", initial_error)
    print("Initial gradient: ", initial_gradient)