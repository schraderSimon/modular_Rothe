"""
Rothe solver for the time-dependent Schrödinger equation.

This module contains the main RotheSolver class for implicit Rothe time stepping
in a non-orthogonal Gaussian basis, plus the factory function for setting up
the Rothe error and gradient computation.
"""

import os
import sys
import time

os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true --xla_gpu_autotune_level=4"

import numpy as np
from scipy.optimize import minimize

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from jax.scipy.linalg import solve

from .file_handling import OutputConfig, _select_and_trim_to_time, save_rothe_state
from .general_wf import generalPotentialSolver, propagate_kinetic_analytical

# Import optimization helpers
from .optimization_helpers import (
    build_diagonal_scaling_from_grad,
    compute_gtol,
    make_initial_hessian_y,
    make_objective_params,
    make_objective_y,
    maybe_apply_bb1_scaling,
    params_flat_to_y_flat,
    y_flat_to_params_flat,
)

# Re-export OutputConfig and save_rothe_state for backward compatibility
__all__ = ["RotheSolver", "setUpRotheErrorAndGradient_jit", "save_rothe_state", "OutputConfig"]


class RotheSolver:
    """
    You provide:
        - `SHH2(t, params, splitting_type) -> (S, H, H2)` where
            S, H, H2 are (N,N) matrices in the Gaussian basis.
        - `rothe_grad_jax`: a JAX function returning (value, grad) w.r.t. p_new.
        - `rothe_nograd`: a non-grad version (used to compute c_new).

    State:
        - `params_old`: Gaussian parameters (typically shape (n, D, 4) with (a,b,mu,p)).
        - `coeffs_old`: coefficients (shape (n,), complex).

    Output:
        - Optionally appends to an HDF5 log via `file_handling.save_rothe_state`.
    """

    def __init__(
        self,
        SHH2,
        dt,
        t,
        params_old,
        coeffs_old,
        rothe_grad_jax,
        rothe_nograd=None,
        splitting_type="none",
        output_config=None,
    ):
        # Core solver state
        self.SHH2 = SHH2
        self.dt = dt
        self.t = t
        self.params_old = jnp.asarray(params_old)
        self.coeffs_old = jnp.asarray(coeffs_old)
        self.rothe_error_and_gradient = rothe_grad_jax
        self.rothe_nograd = rothe_nograd
        self.splitting_type = splitting_type
        self.total_time = 0
        self.params_oldold = None

        # Output configuration
        if output_config is None:
            output_config = OutputConfig()
        self.output_config = output_config
        self.name = output_config.name
        self.polynomial_string = output_config.polynomial_string
        self.out_dir = output_config.out_dir
        self.compression = output_config.compression
        self.compression_opts = output_config.compression_opts

    def find_next_timestep_solution(self, params_init=None, maxiter=100):
        # 1) Choose initial p and store shape
        if params_init is None:
            params_init = np.asarray(self.params_old)
        self.params_init = np.asarray(params_init)

        params_shape = self.params_init.shape
        x0_params_flat = self.params_init.ravel()  # initial point in params-space
        N = x0_params_flat.size

        # 2) p-space objective + gradient (needed in any case)
        f_and_g_params, g_only_params = make_objective_params(
            self.rothe_error_and_gradient,
            params_shape,
            self.params_old,
            self.coeffs_old,
            self.SHH2,
            self.t,
            self.dt,
        )

        # 3) Evaluate at initial guess in p-space
        initial_err, initial_grad_params_flat = f_and_g_params(x0_params_flat)
        print("Initial squared Rothe error: ", initial_err)
        gtol = compute_gtol(initial_err)

        grad_params = initial_grad_params_flat.reshape(params_shape)
        S_flat, S_flat_jax = build_diagonal_scaling_from_grad(grad_params)
        # If you actually want scaling, comment out these two lines:
        # S_flat_jax = jnp.ones_like(S_flat_jax)
        # S_flat = np.ones_like(S_flat)

        # y0 and gradient in y-space
        theta0 = params_flat_to_y_flat(x0_params_flat, S_flat)
        grad_theta0_flat = initial_grad_params_flat * S_flat  # g_y = S g_params
        grad_theta0_norm = np.linalg.norm(grad_theta0_flat)

        # linear map p = T theta with theta=y
        T = np.diag(S_flat)  # p = S y
        T_inv = np.diag(1.0 / S_flat)  # y = S^{-1} p

        # y-space objective
        f_and_g_theta = make_objective_y(
            self.rothe_error_and_gradient,
            params_shape,
            self.params_old,
            self.coeffs_old,
            self.SHH2,
            self.t,
            self.dt,
            S_flat_jax,
        )

        # 4) Initial inverse Hessian (possibly BB1-scaled) in theta-space
        hess_inv_default = make_initial_hessian_y(N, grad_theta0_norm)

        hess_inv0 = maybe_apply_bb1_scaling(
            params_old=self.params_old,
            params_oldold=getattr(self, "params_oldold", None),
            T=T,
            T_inv=T_inv,
            grad_theta0_flat=grad_theta0_flat,
            g_only_params=g_only_params,
            t=self.t,
            hess_inv_default=hess_inv_default,
        )
        # hess_inv0 = np.zeros_like(hess_inv0)
        # print("Gradient norm is: %.3e" % grad_theta0_norm)
        # for i in range(len(hess_inv0)):
        #    hess_inv0[i, i] = 1.0 / grad_theta0_norm
        if maxiter == 0:
            niter = -1
            theta_solved = theta0
        else:
            solution = minimize(
                f_and_g_theta,
                theta0,
                jac=True,
                method="BFGS",
                options={
                    "gtol": gtol,
                    "maxiter": maxiter,
                    "disp": False,
                    "hess_inv0": hess_inv0,
                },
            )
            niter = solution.nit
            theta_solved = solution.x

        params_solved_flat = y_flat_to_params_flat(theta_solved, S_flat)

        params_solved = params_solved_flat.reshape(params_shape)

        # 7) Final Rothe error and new coefficients
        final_RE, coeffs_new = self.rothe_nograd(
            jnp.asarray(params_solved),
            self.params_old,
            self.coeffs_old,
            self.SHH2,
            self.t,
            self.dt,
            return_cnew=True,
        )
        return params_solved, coeffs_new, float(final_RE), niter

    def propagate(
        self, num_iterations, params_init=None, coeffs_init=None, maxiter=100, renormalize=True
    ):
        if params_init is not None:
            self.params_old = jnp.asarray(params_init)
        if coeffs_init is not None:
            self.coeffs_old = jnp.asarray(coeffs_init)
        if self.params_oldold is not None:
            startguess = 2 * self.params_old - self.params_oldold
        else:
            startguess = self.params_old
        D = int(self.params_old.shape[1])
        onesD = jnp.ones(D)
        for i in range(num_iterations):
            self.t += self.dt
            dt = self.dt
            start = time.time()
            if self.splitting_type in (None, "none"):
                params_solved, coeffs_new, final_RE, nit = self.find_next_timestep_solution(
                    startguess, maxiter=maxiter
                )
            elif self.splitting_type == "kinetic":
                # Half step with T

                params_old_temp = self.params_old.copy()
                coeffs_old_temp = self.coeffs_old.copy()
                self.params_old, self.coeffs_old = propagate_kinetic_analytical(
                    dt / 2, self.params_old, self.coeffs_old, onesD
                )
                # Need to update coeffs/params_old as find_next_timestep_solution uses it as its base point
                startguess, _ = propagate_kinetic_analytical(
                    -dt / 2, startguess, self.coeffs_old, onesD
                )
                params_solved_temp, coeffs_new_temp, final_RE, nit = (
                    self.find_next_timestep_solution(startguess, maxiter=maxiter)
                )
                params_solved, coeffs_new = propagate_kinetic_analytical(
                    dt / 2, params_solved_temp, coeffs_new_temp, onesD
                )
                self.params_old = params_old_temp  # Restore to original values
                self.coeffs_old = coeffs_old_temp
            else:
                raise ValueError(
                    f"Unknown splitting_type {self.splitting_type!r}; expected 'none' or 'kinetic'."
                )
            if np.isrealobj(self.t) and not np.iscomplexobj(self.t):
                time_str = f"{float(np.real(self.t)):.2f}"
            else:
                time_str = str(self.t)
            print(
                f"Time {time_str}:Squared Rothe error: {final_RE:.3e}, Number of iterations: {nit}"
            )
            if renormalize:
                S, H, H2 = self.SHH2(
                    t=self.t, params=params_solved, splitting_type=self.splitting_type
                )
                norm = jnp.conj(coeffs_new) @ S @ coeffs_new

                sqrt_norm = jnp.sqrt(norm)
                coeffs_new = coeffs_new / sqrt_norm  # enforce <c|S|c>=1
                if self.dt.imag != 0:
                    print("Norm before renormalization: ", norm)
                    new_energy = jnp.real(jnp.conj(coeffs_new) @ H @ coeffs_new)
                    H2_expec = jnp.real(jnp.conj(coeffs_new) @ H2 @ coeffs_new)

                    print("Energy after renormalization: ", new_energy)
                    variance = H2_expec - new_energy**2
                    print("Variance after renormalization: ", variance)
            if self.name is not None:
                # Calculate dipole moment if polynomial_string is available
                dipole_moment = None
                if self.polynomial_string is not None:
                    try:
                        solver = generalPotentialSolver(params_solved, self.polynomial_string)
                        dipole_moment = solver.calculate_dipole_moment(coeffs_new)
                    except Exception as e:
                        print(f"Warning: Could not compute dipole moment: {e}")

                save_rothe_state(
                    self.name,
                    self.splitting_type,
                    self.polynomial_string or "",
                    self.dt,
                    self.t,
                    final_RE,
                    params_solved,
                    coeffs_new,
                    dipole_moment=dipole_moment,
                    compression=self.compression,
                    compression_opts=self.compression_opts,
                    path=self.out_dir,
                )
            if self.params_old is None:
                startguess = params_solved
            startguess = 2 * params_solved - self.params_old

            self.params_oldold = self.params_old
            self.params_old = params_solved

            self.coeffs_old = coeffs_new
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
        if sel["params"].shape != tuple(self.params_old.shape):
            raise ValueError(
                f"Saved params shape {sel['params'].shape} != current {tuple(self.params_old.shape)}"
            )
        if sel["coeffs"].shape != tuple(self.coeffs_old.shape):
            raise ValueError(
                f"Saved coeffs shape {sel['coeffs'].shape} != current {tuple(self.coeffs_old.shape)}"
            )

        self.params_old = jnp.asarray(sel["params"])
        self.coeffs_old = jnp.asarray(sel["coeffs"])
        self.t = np.asarray(sel["t"]).item()
        if sel["params_prev"] is not None:
            self.params_oldold = jnp.asarray(sel["params_prev"])
        else:
            self.params_oldold = None
        return sel  # contains t, idx, trimmed flag


def setUpRotheErrorAndGradient_jit(splitting_type):
    def rothe_error(
        params_new, params_old, coeffs_old, SHH2, t, dt, return_cnew=False, lambda_=1e-10
    ):
        splitting_eff = "none" if splitting_type in (None, "none") else splitting_type
        ngo = params_old.shape[0]
        params_concat = jnp.concatenate(
            (params_old, params_new), axis=0
        )  # New shape: (2*n,N_b)
        S_full, H_full, H2_full = SHH2(
            t=t + dt / 2, params=params_concat, splitting_type=splitting_eff
        )
        dtc = jnp.conj(dt)
        dt_abs_sq_4 = dt * dtc / 4

        A_dagger_A = (
            S_full[ngo:, ngo:]
            + dt_abs_sq_4 * H2_full[ngo:, ngo:]
            + 0.5j * (dt - dtc) * H_full[ngo:, ngo:]
        )

        B_dagger_B = (
            S_full[:ngo, :ngo]
            + dt_abs_sq_4 * H2_full[:ngo, :ngo]
            - 0.5j * (dt - dtc) * H_full[:ngo, :ngo]
        )

        rho_mat = (
            S_full[:ngo, ngo:]
            - dt_abs_sq_4 * H2_full[:ngo, ngo:]
            - 0.5j * (dt + dtc) * H_full[:ngo, ngo:]
        )
        rho_vec = jnp.conj(rho_mat).T @ coeffs_old
        S_reg = A_dagger_A + jnp.eye(A_dagger_A.shape[0]) * lambda_
        coeffs_new = solve(S_reg, rho_vec)

        overlap_term = jnp.conj(coeffs_old) @ B_dagger_B @ coeffs_old
        projection_term = jnp.conj(rho_vec).T @ coeffs_new
        rothe_error = jnp.real(overlap_term - projection_term)
        if return_cnew:
            return rothe_error, coeffs_new
        else:
            return rothe_error

    rothe_vg_jit = jax.jit(
        jax.value_and_grad(rothe_error, argnums=0), static_argnames=("SHH2",)
    )
    return rothe_error, rothe_vg_jit
