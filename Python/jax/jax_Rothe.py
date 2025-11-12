import os
import time

import h5py
import numpy as np

# set these BEFORE importing jax
from scipy.optimize import minimize

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)
import jax.scipy as jsp
from file_handling import _select_and_trim_to_time, load_rothe_file, save_rothe_state
from general_wf import propagate_kinetic_analytical
from jax.scipy.linalg import solve


class RotheSolver:
    """
    This class implements the Rothe method for the time-dependent Schrödinger equation.
    Essentially, it does the following:
    1. The user has to provide a function f()
    """

    def __init__(
        self,
        SHH2,
        dt,
        t,
        p_old,
        c_old,
        rothe_grad_jax,
        rothe_nograd=None,
        splitting_type="none",
        name=None,
        polynomial_string=None,
        out_dir=".",
        compression="gzip",
        compression_opts=4,
    ):
        self.SHH2 = SHH2
        self.dt = dt
        self.t = t
        self.p_old = jnp.asarray(p_old)
        self.c_old = jnp.asarray(c_old)
        self.rothe_error_and_gradient = rothe_grad_jax
        self.rothe_nograd = rothe_nograd
        self.splitting_type = splitting_type
        self.total_time = 0
        self.p_oldold = None
        self.name = name
        self.polynomial_string = polynomial_string
        self.out_dir = out_dir
        self.compression = compression
        self.compression_opts = compression_opts

    def find_next_timestep_solution(self, p_init=None, maxiter=100):
        if p_init is None:
            p_init = np.asarray(self.p_old)
        self.p_init = np.asarray(p_init)

        def f_and_g(theta_flat):
            p = jnp.asarray(theta_flat).reshape(self.p_init.shape)
            val, g = self.rothe_error_and_gradient(
                p, self.p_old, self.c_old, self.SHH2, self.t, self.dt
            )
            return float(val), np.asarray(jnp.ravel(g))

        def f_only(theta):
            v, _ = f_and_g(theta)
            return v

        def g_only(theta):
            _, g = f_and_g(theta)
            return g

        x0 = self.p_init.ravel()
        # hess_inv=np.eye(len(x0)) #Initial Hessian
        initial_err, initial_gradient = f_and_g(x0)
        print("Initial Rothe error: ", initial_err)
        hess_inv = np.eye(len(x0)) / jnp.linalg.norm(initial_gradient)

        if self.p_oldold is not None:
            gradient_oldold = g_only(jnp.ravel(self.p_oldold))
            y = initial_gradient - gradient_oldold
            s = jnp.ravel(self.p_old) - jnp.ravel(self.p_oldold)
            BB1 = jnp.inner(s, s) / (jnp.inner(s, y))
            if BB1 > 0:  # Make sure curvature condition is satisfied
                print("%f fulfills curvature condition: BB1=%f" % (self.t, BB1))
                hess_inv = BB1 * jnp.eye(len(x0))
            else:
                print("%f does not fulfill curvature condition: BB1=%f" % (self.t, BB1))
        solution = minimize(
            f_and_g,
            x0,
            jac=True,
            method="BFGS",
            options={"gtol": 1e-12, "maxiter": maxiter, "disp": False, "hess_inv0": hess_inv},
        )
        niter = solution.nfev
        p_solved = solution.x.reshape(self.p_init.shape)

        final_RE, c_new = self.rothe_nograd(
            jnp.asarray(p_solved),
            self.p_old,
            self.c_old,
            self.SHH2,
            self.t,
            self.dt,
            return_cnew=True,
        )
        return p_solved, c_new, float(final_RE), niter

    """
    def find_next_timestep_solution(self, p_init=None, maxiter=100, history_size=10, tol=1e-12):
        if p_init is None:
            p_init = self.p_old
        p0 = jnp.asarray(p_init)

        # rothe_error_and_gradient already returns (value, grad) w.r.t. p_new
        def fun(p_new):
            val, g = self.rothe_error_and_gradient(
                p_new, self.p_old, self.c_old, self.SHH2, self.t, self.dt
            )
            return val, g  # grad has same shape as p_new

        # L-BFGS with strong-Wolfe line search (“zoom”) and correct gamma-scaling
        solver = LBFGS(
            fun=fun,
            value_and_grad=True,
            maxiter=maxiter,
            tol=tol,
            history_size=history_size,
            use_gamma=True,        # H0 = γ I with γ from Nocedal & Wright (7.20)
        )

        params, state = solver.run(init_params=p0)   # fully JIT’d on GPU
        p_solved = params

        final_RE, c_new = self.rothe_nograd(
            p_solved, self.p_old, self.c_old, self.SHH2, self.t, self.dt, return_cnew=True
        )
        niter = int(state.iter_num)
        return np.asarray(p_solved), np.asarray(c_new), float(final_RE), niter
    """

    def propagate(self, num_iterations, p_init=None, c_init=None, maxiter=100):
        if p_init is not None:
            self.p_old = jnp.asarray(p_init)
        if c_init is not None:
            self.c_old = jnp.asarray(c_init)
        startguess = self.p_old
        D = int(self.p_old.shape[1])
        onesD = jnp.ones(D)
        for i in range(num_iterations):
            self.t += self.dt
            dt = self.dt
            start = time.time()
            if self.splitting_type == "none":  # If no splitting is used
                p_solved, c_new, final_RE, nit = self.find_next_timestep_solution(
                    startguess, maxiter=maxiter
                )
            elif (
                self.splitting_type == "kinetic"
            ):  # If analytical propagation and splitting is used
                # Half step with T
                self.p_old, self.c_old = propagate_kinetic_analytical(
                    dt / 2, self.p_old, self.c_old, onesD
                )
                startguess, _ = propagate_kinetic_analytical(
                    -dt / 2, startguess, self.c_old, onesD
                )
                p_solved_temp, c_new_temp, final_RE, nit = self.find_next_timestep_solution(
                    startguess, maxiter=maxiter
                )
                p_solved, c_new = propagate_kinetic_analytical(
                    dt / 2, p_solved_temp, c_new_temp, onesD
                )
            print(
                f"Time {self.t:.2f}: Rothe error: {final_RE:.3e}, Number of iterations: {nit}"
            )
            if self.name is not None:
                save_rothe_state(
                    self.name,
                    self.splitting_type,
                    self.polynomial_string or "",
                    self.dt,
                    self.t,
                    final_RE,
                    p_solved,
                    c_new,
                    compression=self.compression,
                    compression_opts=self.compression_opts,
                    path=self.out_dir,
                )
            startguess = 2 * p_solved - self.p_old
            startguess = p_solved
            self.p_oldold = self.p_old
            self.p_old = p_solved
            self.c_old = c_new
            end = time.time()
            time_taken = end - start
            self.total_time += time_taken

            print(f"Time taken: {time_taken:.2f}")
            # print(p_solved)
            # print(c_new)
        print(f"Total time: {self.total_time}")

    def resume_from_file(self, t_start=None):
        """
        Load state from HDF5 and set solver to resume at requested time.
        If t_start is None -> resume from last saved step.
        If t_start exists -> resume at that step and delete later steps.
        If t_start doesn't exist -> resume from latest step < t_start and delete all >= t_start.
        """
        if self.name is None:
            raise ValueError("Cannot resume: 'name' is None (no filename to load).")
        filename = os.path.join(self.out_dir, f"{self.name}__{self.splitting_type}.h5")
        sel = _select_and_trim_to_time(filename, t_start)

        # optional safety: check shapes
        if sel["p"].shape != tuple(self.p_old.shape):
            raise ValueError(
                f"Saved p shape {sel['p'].shape} != current {tuple(self.p_old.shape)}"
            )
        if sel["c"].shape != tuple(self.c_old.shape):
            raise ValueError(
                f"Saved c shape {sel['c'].shape} != current {tuple(self.c_old.shape)}"
            )

        self.p_old = jnp.asarray(sel["p"])
        self.c_old = jnp.asarray(sel["c"])
        self.t = float(sel["t"])
        self.p_oldold = None  # reset history across resume
        return sel  # contains t, idx, trimmed flag


def setUpRotheErrorAndGradient_jit(splitting_type):
    def rothe_error(p_new, p_old, c_old, SHH2, t, dt, return_cnew=False, lambda_=1e-8):
        ngo = p_old.shape[0]
        dtsq4 = dt**2 / 4
        p_concat = jnp.concatenate((p_old, p_new), axis=0)  # New shape: (2*n,N_b)
        S_full, H_full, H2_full = SHH2(
            t=t + dt / 2, params=p_concat, splitting_type=splitting_type
        )
        S_tilde_full = S_full + dtsq4 * H2_full
        rho_mat = (
            S_full[:ngo, ngo:] - dtsq4 * H2_full[:ngo, ngo:] + 1j * dt * H_full[:ngo, ngo:]
        )
        rho_vec = jnp.conj(rho_mat).T @ c_old
        S_tilde_full = 0.5 * (S_tilde_full + jnp.conj(S_tilde_full.T))
        S_reg = (
            S_tilde_full + jnp.eye(S_tilde_full.shape[0]) * lambda_
        )  # Regularization term to avoid singular matrices.
        # L = jnp.linalg.cholesky(S_reg[ngo:, ngo:])  # Hermitian PD
        # c_new = jsp.linalg.solve_triangular(L, rho_vec, lower=True, trans="N")
        # c_new = jsp.linalg.solve_triangular(L.conj().T, c_new, lower=False, trans="N")
        c_new = solve(S_reg[ngo:, ngo:], rho_vec, assume_a="her")
        overlap_term = jnp.conj(c_old) @ S_tilde_full[:ngo, :ngo] @ c_old
        projection_term = jnp.real(jnp.conj(rho_vec).T @ c_new)
        rothe_error = jnp.real(overlap_term - projection_term)
        if return_cnew:
            return rothe_error, c_new
        else:
            return rothe_error

    rothe_vg_jit = jax.jit(
        jax.value_and_grad(rothe_error, argnums=0), static_argnames=("SHH2",)
    )
    return rothe_error, rothe_vg_jit
