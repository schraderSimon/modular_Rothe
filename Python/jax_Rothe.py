import matplotlib.pyplot as plt
import jax.numpy as jnp
import jax
from jax.scipy.linalg import solve
from scipy.optimize import minimize
import h5py
import numpy as np
import sys
import os
import time
# set these BEFORE importing jax
from jaxopt import BFGS
jax.config.update('jax_enable_x64', True)

def rothe_error(p_new,p_old,c_old,SHH2,t,dt,return_cnew=False,lambda_=1e-14):

    ngo= p_old.shape[0]
    ngn= p_new.shape[0]

    p_concat=jnp.concatenate((p_old,p_new),axis=0) #New shape: (2*n,N_b)
    S_full,H_full,H2_full=SHH2(t=t+dt/2,params=p_concat)
    S_tilde_full=S_full + dt**2/4* H2_full
    rho_mat   = (S_full[:ngo, ngo:] - (dt**2/4) * H2_full[:ngo, ngo:] + 1j*dt*H_full[:ngo, ngo:])
    rho_vec   = jnp.conj(rho_mat).T @ c_old
    S_tilde_full=0.5*(S_tilde_full + jnp.conj(S_tilde_full.T))
    S_reg= S_tilde_full + jnp.eye(S_tilde_full.shape[0])*lambda_ #Regularization term to avoid singular matrices.
    c_new = solve(S_reg[ngo:,ngo:],rho_vec,assume_a="her")
    overlap_term = jnp.conj(c_old)@S_tilde_full[:ngo,:ngo]@c_old
    projection_term = jnp.real(jnp.conj(rho_vec).T @ c_new )
    rothe_error=jnp.real(overlap_term - projection_term)
    if return_cnew:
        return rothe_error, c_new
    else:
        return rothe_error
rothe_vg_jit = jax.jit(jax.value_and_grad(rothe_error, argnums=0), static_argnames=('SHH2',))
class RotheSolver:
    """
    This class implements the Rothe method for the time-dependent Schrödinger equation.
    Essentially, it does the following:
    1. The user has to provide a function f()
    """
    def __init__(self, SHH2, dt,t,p_old,c_old):
        self.SHH2 = SHH2
        self.dt = dt
        self.t = t
        self.p_old = jnp.asarray(p_old)
        self.c_old = jnp.asarray(c_old)
    def find_next_timestep_solution(self, p_init=None,maxiter=100):
        if p_init is None:
            p_init = np.asarray(self.p_old)
        self.p_init = np.asarray(p_init)

        def f_and_g(theta_flat):
            p = jnp.asarray(theta_flat).reshape(self.p_init.shape)
            val, g = rothe_vg_jit(p, self.p_old, self.c_old, self.SHH2, self.t, self.dt)
            return float(val), np.asarray(jnp.ravel(g))

        def f_only(theta): v, _ = f_and_g(theta); return v
        def g_only(theta): _, g = f_and_g(theta); return g

        x0 = self.p_init.ravel()
        hess_inv=np.eye(len(x0)) #Initial Hessian
        #hess_inv=np.diag(1/(abs(g_only(x0)+1e-14))) #Preconditioner
        
        solution = minimize(f_only, x0, jac=g_only, method='BFGS',
                            options={'gtol': 1e-14, 'maxiter': maxiter,"disp":True,"hess_inv0":hess_inv})
        niter = solution.nfev
        p_solved = solution.x.reshape(self.p_init.shape)

        final_RE, c_new = rothe_error(
            jnp.asarray(p_solved), self.p_old, self.c_old, self.SHH2, self.t, self.dt, return_cnew=True
        )
        return p_solved, c_new, float(final_RE), niter
    def propagate(self,num_iterations,p_init=None,c_init=None):
        if p_init is not None:
            self.p_old = jnp.asarray(p_init)
        if c_init is not None:
            self.c_old = jnp.asarray(c_init)
        startguess = self.p_old
        for i in range(num_iterations):
            self.t += self.dt
            start=time.time()
            p_solved, c_new, final_RE, nit = self.find_next_timestep_solution(startguess)
            print(f"Step {i}: Rothe error: {final_RE}, Number of iterations: {nit}")
            startguess = 2*p_solved-self.p_old
            self.p_old = p_solved
            self.c_old = c_new
            end=time.time()
            print(f"Time taken for step {i}: {end-start}")
if __name__ == "__main__":
    from HOscillator_jax import HOscillator_ND
    # ---------- initial wave-function parameters ------------------------
    dt     = 0.01
    omega  = 1.0
    n=3
    D=2
    def make_SHH2_pure(omega=1.0, K_max=4):
        def SHH2(t=0.0, params=None):
            osc = HOscillator_ND(params, omega=omega, K_max=K_max)  # new, clean object every time
            return osc.calculate_SHH2(t=t)
        return SHH2
    p_init = jnp.zeros((n, D, 4),dtype=jnp.float64)
    p_init = p_init.at[:, :, 0].set(1/jnp.sqrt(2)) #Width parameters
    for d in range(D):
        p_init = p_init.at[:, d, 2].set(jnp.linspace(-2,2,n,dtype=jnp.float64))          #Move to the right
        #p_init = p_init.at[:, d, 3].set(jnp.linspace(2,-2,n,dtype=jnp.float64))          #Set p to 1
    # ---------- build WF + solver ---------------------------------------
    oscillator=HOscillator_ND(p_init,omega=omega)
    SHH2=make_SHH2_pure(omega) # The function to calculate S, H, H2
    Smat = oscillator.calculate_S()
    c_init = jnp.array([1.0]*n,dtype=jnp.complex128)  # c
    norm=c_init.T@Smat@c_init
    c_init=c_init/jnp.sqrt(norm)
    initial_error=rothe_error(p_init,p_init,c_init,SHH2,0,dt)
    gradient_function=grad_err = jax.grad(rothe_error, argnums=0)
    initial_gradient=gradient_function(p_init,p_init.copy(),c_init,SHH2,0,dt)
    print("Initial Rothe error: ", initial_error)
    print("Initial Rothe gradient: ", initial_gradient)
    rothesolver=RotheSolver(SHH2, dt, 0, p_init, c_init)
    p_solved,c_new,final_RE,nit=rothesolver.find_next_timestep_solution(maxiter=1)
    print("Rothe error: ", final_RE)
    print("Number of iterations: ", nit)
    
    p_init=p_init+np.random.uniform(-0.01,0.01,p_init.shape)
    start=time.time()
    rothesolver=RotheSolver(SHH2, dt, 0, p_init, c_init)
    print("Initial Rothe error: %e"%(rothe_error(p_init,p_init,c_init,SHH2,0,dt)))
    p_solved,c_new,final_RE,nit=rothesolver.find_next_timestep_solution(maxiter=30)
    print("Rothe error: ", final_RE)
    print("Number of iterations: ", nit)
    end=time.time()
    print("Time taken: ", end-start)
    print("Time per iteration: ", (end-start)/nit if nit > 0 else 0)
    rothesolver.propagate(200,p_init,c_init)