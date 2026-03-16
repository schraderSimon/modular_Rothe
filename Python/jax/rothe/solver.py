"""
Rothe solver for the time-dependent Schrödinger equation.

This module contains the main RotheSolver class for implicit Rothe time stepping
in a non-orthogonal Gaussian basis, plus the factory function for setting up
the Rothe error and gradient computation.
"""

import os
import time
from dataclasses import dataclass
from functools import partial

import numpy as np
from scipy.optimize import minimize

import jax
import jax.numpy as jnp
from rothe.io import OutputConfig, _select_and_trim_to_time, load_rothe_file, save_rothe_state
from rothe.objective import (
    _log_det_penalty,
    _overlap_penalty,
    build_objective_context,
    compute_SHH2_oldold,
    setUpRotheErrorAndGradient_jit,
)
from rothe.optimization import (
    compute_gtol,
    compute_tanh_bounds,
    make_objective_params,
    make_objective_tanh,
    maybe_apply_bb1_scaling,
    tanh_forward,
    tanh_inverse,
)
from rothe.wavefunction import generalPotentialSolver, propagate_kinetic_analytical

__all__ = ["RotheSolver", "setUpRotheErrorAndGradient_jit", "save_rothe_state", "OutputConfig"]


_CANDIDATE_OVERLAP_DISCARD_THRESHOLD = 0.99


@partial(jax.jit, static_argnames=("SHH2", "splitting_eff", "has_frozen"))
def _augment_blocks_with_candidate(
    SHH2,
    t_mid,
    dt_real,
    dt_imag,
    dt_abs_sq_4,
    splitting_eff,
    A_base,
    rho_base,
    params_new_base,
    params_old,
    candidate_param,
    has_frozen,
    params_frozen=None,
    A_df=None,
    A_ff=None,
    rho_xf=None,
    S_base=None,
    S_df=None,
):
    """Assemble augmented blocks for one candidate Gaussian (JIT-compiled).

    Returns (A_aug, rho_aug, S_dyn_aug, S_of_aug) where S_dyn_aug is the
    new dynamic-dynamic overlap matrix (including the candidate) and S_of_aug
    is the new dynamic-frozen overlap (shape (0,0) when has_frozen is False).
    """
    S_cc, H_cc, H2_cc = SHH2(
        t=t_mid, params_ket=candidate_param, params_bra=candidate_param, splitting_type=splitting_eff
    )
    A_cc = S_cc + dt_abs_sq_4 * H2_cc - dt_imag * H_cc

    S_cd, H_cd, H2_cd = SHH2(
        t=t_mid, params_ket=params_new_base, params_bra=candidate_param, splitting_type=splitting_eff
    )
    A_cd = S_cd + dt_abs_sq_4 * H2_cd - dt_imag * H_cd

    S_dc, H_dc, H2_dc = SHH2(
        t=t_mid, params_ket=candidate_param, params_bra=params_new_base, splitting_type=splitting_eff
    )
    A_dc = S_dc + dt_abs_sq_4 * H2_dc - dt_imag * H_dc

    S_xc, H_xc, H2_xc = SHH2(t=t_mid, params_ket=candidate_param, params_bra=params_old, splitting_type=splitting_eff)
    rho_xc = S_xc - dt_abs_sq_4 * H2_xc - 1j * dt_real * H_xc

    if not has_frozen:
        A_aug = jnp.block([[A_base, A_dc], [A_cd, A_cc]])
        rho_aug = jnp.concatenate([rho_base, rho_xc], axis=1)
        S_dyn_aug = jnp.block([[S_base, S_dc], [S_cd, S_cc]])
        S_of_aug = jnp.zeros((0, 0), dtype=S_dyn_aug.dtype)
    else:
        n_dyn = params_new_base.shape[0]
        A_dd = A_base[:n_dyn, :n_dyn]
        rho_xd = rho_base[:, :n_dyn]

        S_cf, H_cf, H2_cf = SHH2(
            t=t_mid, params_ket=params_frozen, params_bra=candidate_param, splitting_type=splitting_eff
        )
        A_cf = S_cf + dt_abs_sq_4 * H2_cf - dt_imag * H_cf

        S_fc, H_fc, H2_fc = SHH2(
            t=t_mid, params_ket=candidate_param, params_bra=params_frozen, splitting_type=splitting_eff
        )
        A_fc = S_fc + dt_abs_sq_4 * H2_fc - dt_imag * H_fc

        A_aug = jnp.block([[A_dd, A_dc, A_df], [A_cd, A_cc, A_cf], [jnp.conj(A_df).T, A_fc, A_ff]])
        rho_aug = jnp.concatenate([rho_xd, rho_xc, rho_xf], axis=1)
        S_dyn_aug = jnp.block([[S_base, S_dc], [S_cd, S_cc]])
        S_of_aug = jnp.concatenate([S_df, S_cf], axis=0)  # shape (n_dyn+1, n_frozen)

    return A_aug, rho_aug, S_dyn_aug, S_of_aug


def _flip_bra_phase_params(params):
    """Return params with bra-phase signs flipped (b -> -b, p -> -p)."""
    params_flipped = jnp.array(params, copy=True)
    params_flipped = params_flipped.at[:, :, 1].set(-params_flipped[:, :, 1])
    params_flipped = params_flipped.at[:, :, 3].set(-params_flipped[:, :, 3])
    return params_flipped


@dataclass(frozen=True)
class StepContext:
    """Immutable context for a single Rothe time step."""

    n_frozen: int
    n_dynamic: int
    ones_D: jnp.ndarray


@dataclass(frozen=True)
class StepResult:
    """Outcome of solving a single Rothe time step."""

    params_solved: jnp.ndarray
    coeffs_new: jnp.ndarray
    rothe_error: float
    n_iterations: int


@dataclass(frozen=True)
class StepMetrics:
    """Observables computed after a time step."""

    max_overlap: float
    dipole_moment: jnp.ndarray | None = None
    double_autocorrelation: jnp.ndarray | None = None


class RotheSolver:
    """
    Rothe solver for the TDSE in a Gaussian basis.

    You provide:
        - ``SHH2(t, params, splitting_type) -> (S, H, H2)`` where
          S, H, H2 are (N,N) matrices in the Gaussian basis.
                - ``rothe_grad_fn``: function returning (value, grad) w.r.t. p_new.
                    It can be JAX-based, NumPy-based, or any callable with matching signature.
        - ``rothe_nograd``: a non-grad version (used to compute c_new).

    State:
        - ``params_old``: **full** Gaussian parameters at the current time — i.e.
          ``concat(optimisable, frozen)`` once ``params_frozen`` is set.
          Shape ``(n_old + n_frozen, D, 4)``.
        - ``coeffs_old``: coefficients (shape ``(n_old + n_frozen,)``), complex.
        - ``params_frozen``: optional shape ``(n_frozen, D, 4)`` array of basis
          functions that are appended to every new-time wave function but never
          optimised. The optimisation variable has shape ``(n_old, D, 4)`` only.

    Output:
        - Optionally appends to an HDF5 log via :func:`save_rothe_state`.
    """

    def __init__(
        self,
        SHH2,
        dt,
        t,
        epsilon,
        regularization_lambda,
        params_old,
        coeffs_old,
        rothe_grad_fn,
        rothe_nograd=None,
        splitting_type="none",
        output_config=None,
        params_frozen=None,
        random_seed=0,
        predictive_penalty_multiplier=0.1,
    ):
        self.SHH2 = SHH2
        self.dt = dt
        self.t = t
        self.epsilon = float(epsilon)
        self.regularization_lambda = float(regularization_lambda)
        self.predictive_penalty_multiplier = float(predictive_penalty_multiplier)
        self.random_seed = int(random_seed)
        self._rng = np.random.default_rng(self.random_seed)
        self.threshold = None
        self.params_old = jnp.asarray(params_old)
        self.coeffs_old = jnp.asarray(coeffs_old)
        self.rothe_grad_fn = rothe_grad_fn
        self.rothe_nograd = rothe_nograd
        self.splitting_type = splitting_type
        self.total_time = 0
        self.params_oldold = None
        self.params_frozen = jnp.asarray(params_frozen) if params_frozen is not None else None
        self.cumulative_rothe_error = 0.0
        self.previous_step_rothe_error = None
        self.gaussian_addition_cooldown_steps = 10
        self._gaussian_addition_cooldown_remaining = 0

        if output_config is None:
            output_config = OutputConfig()
        self.output_config = output_config

        if self.output_config.polynomial_string is not None:
            self._dipole_solver = generalPotentialSolver(
                params_ket=params_old, params_bra=params_old, polynomial_string=self.output_config.polynomial_string
            )
        else:
            self._dipole_solver = None

    def _eval_rothe_nograd(self, params_new, return_cnew=False, SHH2_oldold=None):
        if self.rothe_nograd is None:
            raise ValueError("rothe_nograd must be set to evaluate Rothe objective without gradients.")
        return self.rothe_nograd(
            jnp.asarray(params_new),
            self.params_old,
            self.coeffs_old,
            self.SHH2,
            self.t,
            self.dt,
            self.params_frozen,
            return_cnew=return_cnew,
            lambda_=self.regularization_lambda,
            SHH2_oldold=SHH2_oldold,
        )

    def _prepare_timestep_optimization_inputs(self, params_init):
        n_frozen = self.params_frozen.shape[0] if self.params_frozen is not None else 0
        n_dynamic = self.params_old.shape[0] - n_frozen
        params_dynamic = self.params_old[:n_dynamic]

        if params_init is None:
            params_init = np.asarray(params_dynamic)
        params_init = np.asarray(params_init)

        params_shape = params_init.shape
        x0_params_flat = params_init.ravel()
        return n_dynamic, params_shape, x0_params_flat

    def _build_parameter_optimization_problem(
        self, x0_params_flat, params_shape, objective_context, tanh_p=0.5, tanh_q=0.1
    ):
        f_and_g_params, g_only_params = make_objective_params(self.rothe_grad_fn, params_shape, objective_context)

        initial_err, initial_grad_params_flat = f_and_g_params(x0_params_flat)
        if not np.isfinite(initial_err):
            print("Initial squared Rothe error is non-finite; using fallback objective scale for this step.")
            initial_err = 1.0
        print("Initial squared Rothe error: ", initial_err)
        gtol = compute_gtol(initial_err)

        predictive_penalty_scale = 0.0
        if self.previous_step_rothe_error is not None:
            prev_err = float(self.previous_step_rothe_error)
            denom = max(abs(float(initial_err)), 1e-30)
            if np.isfinite(prev_err) and prev_err > 0.0:
                predictive_penalty_scale = self.predictive_penalty_multiplier * prev_err / denom

        penalized_rothe_vg = partial(
            self.rothe_grad_fn,
            overlap_penalty_scale=float(initial_err),
            predictive_penalty_scale=float(predictive_penalty_scale),
        )
        f_and_g_params_penalized, _ = make_objective_params(penalized_rothe_vg, params_shape, objective_context)

        # Tanh-bounded transform: optimise in unconstrained α' space.
        lb, ub = compute_tanh_bounds(x0_params_flat, params_shape, p=tanh_p, q=tanh_q)
        f_and_g_tanh = make_objective_tanh(f_and_g_params_penalized, lb, ub)
        # Initial point in transformed space is always 0.
        x0_prime = np.zeros_like(x0_params_flat)

        N = x0_params_flat.size
        # Scale-invariant default inverse-Hessian seed for BFGS.
        grad_norm = float(np.linalg.norm(initial_grad_params_flat))
        if not np.isfinite(grad_norm) or grad_norm <= 0.0:
            hess_scale = 1.0
        else:
            hess_scale = 1.0 / grad_norm
        hess_inv_default = hess_scale * np.eye(N)
        return (f_and_g_tanh, g_only_params, x0_prime, lb, ub, hess_inv_default, gtol)

    def find_next_timestep_solution(self, params_init=None, maxiter=100, given_gtol=None):
        n_dynamic, params_shape, x0_params_flat = self._prepare_timestep_optimization_inputs(params_init)
        objective_context = build_objective_context(
            self.SHH2,
            self.params_old,
            self.coeffs_old,
            self.t,
            self.dt,
            self.splitting_type,
            self.params_frozen,
            self.regularization_lambda,
        )
        f_and_g_tanh, g_only_params, x0_prime, lb, ub, hess_inv_default, gtol = (
            self._build_parameter_optimization_problem(x0_params_flat, params_shape, objective_context)
        )
        if given_gtol is not None:
            gtol = given_gtol
        params_dynamic = self.params_old[:n_dynamic]
        # Skip BB1 if params_init has more Gaussians than params_old (Gaussian was just added),
        # or if params_oldold has a different dynamic count (Gaussian added at previous step).
        n_init = params_shape[0]
        if self.params_oldold is not None and n_init == n_dynamic and self.params_oldold.shape[0] == n_dynamic:
            params_dynamic_prev = self.params_oldold
        else:
            params_dynamic_prev = None
        hess_inv0 = maybe_apply_bb1_scaling(
            params_old=params_dynamic,
            params_oldold=params_dynamic_prev,
            g_only_params=g_only_params,
            t=self.t,
            hess_inv_default=hess_inv_default,
        )

        if maxiter == 0:
            niter = -1
            alpha_prime_solved = x0_prime
        else:
            solution = minimize(
                f_and_g_tanh,
                x0_prime,
                jac=True,
                method="BFGS",
                options={"gtol": gtol, "maxiter": maxiter, "disp": False, "hess_inv0": hess_inv0},
            )
            niter = solution.nit
            alpha_prime_solved = solution.x

        params_solved_flat = tanh_forward(alpha_prime_solved, lb, ub)
        params_solved = params_solved_flat.reshape(params_shape)

        final_RE, coeffs_new = self._eval_rothe_nograd(
            params_solved, return_cnew=True, SHH2_oldold=objective_context.SHH2_oldold
        )

        # If optimization made things worse, fall back to the start guess.
        params_start = x0_params_flat.reshape(params_shape)
        init_RE, coeffs_start = self._eval_rothe_nograd(
            params_start, return_cnew=True, SHH2_oldold=objective_context.SHH2_oldold
        )
        if float(init_RE) < float(final_RE):
            print(
                f"  Optimization worsened error ({float(final_RE):.3e} > {float(init_RE):.3e}), "
                "keeping initial parameters."
            )
            params_solved = params_start
            coeffs_new = coeffs_start
            final_RE = float(init_RE)

        return params_solved, coeffs_new, float(final_RE), niter

    def init_step_context(self, params_init=None, coeffs_init=None):
        if params_init is not None:
            self.params_old = jnp.asarray(params_init)
        if coeffs_init is not None:
            self.coeffs_old = jnp.asarray(coeffs_init)

        n_frozen = self.params_frozen.shape[0] if self.params_frozen is not None else 0
        n_dynamic = self.params_old.shape[0] - n_frozen
        D = int(self.params_old.shape[1])
        return StepContext(n_frozen=n_frozen, n_dynamic=n_dynamic, ones_D=jnp.ones(D))

    def build_startguess(self, step_context):
        params_dynamic = np.asarray(self.params_old[: step_context.n_dynamic])
        if self.params_oldold is None:
            return params_dynamic

        params_dynamic_prev = np.asarray(self.params_oldold)
        n_common = min(params_dynamic_prev.shape[0], step_context.n_dynamic)
        if n_common <= 0:
            return params_dynamic

        # Build oldold candidate: use params_oldold for first n_common, params_old for the rest.
        startguess_oldold = np.array(params_dynamic, copy=True)
        startguess_oldold[:n_common] = params_dynamic_prev[:n_common]

        # Build extrapolated candidate: 2*old - oldold for first n_common, params_old for the rest.
        startguess_extrapolated = np.array(params_dynamic, copy=True)
        startguess_extrapolated[:n_common] = 2 * params_dynamic[:n_common] - params_dynamic_prev[:n_common]
        # Width parameter a must stay positive.
        startguess_extrapolated[:, :, 0] = np.maximum(np.abs(startguess_extrapolated[:, :, 0]), 1e-6)

        SHH2_oldold = compute_SHH2_oldold(self.SHH2, self.params_old, self.t + self.dt / 2, self.splitting_type)

        candidates = [
            ("params_old", params_dynamic),
            ("params_oldold", startguess_oldold),
            ("params_extrapolated", startguess_extrapolated),
        ]
        errors = [float(self._eval_rothe_nograd(c, SHH2_oldold=SHH2_oldold)) for _, c in candidates]

        finite_indices = [i for i, e in enumerate(errors) if np.isfinite(e)]
        if not finite_indices:
            print("  All start guess objectives non-finite, falling back to params_old.")
            return params_dynamic

        best_idx = min(finite_indices, key=lambda i: errors[i])
        best_name, best_guess = candidates[best_idx]

        other_errs = ", ".join(f"{candidates[i][0]}: {errors[i]:.3e}" for i in finite_indices if i != best_idx)
        if other_errs:
            print(f"  Best start guess: {best_name} ({errors[best_idx]:.3e}), others: {other_errs}.")
        else:
            print(f"  Best start guess: {best_name} ({errors[best_idx]:.3e}).")

        return best_guess

    def solve_step_none(self, startguess, maxiter=100):
        params_solved, coeffs_new, final_RE, nit = self.find_next_timestep_solution(startguess, maxiter=maxiter)
        return StepResult(
            params_solved=params_solved, coeffs_new=coeffs_new, rothe_error=float(final_RE), n_iterations=nit
        )

    def solve_step_kinetic(self, start_guess, step_context, maxiter=100):
        dt = self.dt
        params_old_temp = self.params_old.copy()
        coeffs_old_temp = self.coeffs_old.copy()
        params_frozen_next = self.params_frozen

        self.params_old, self.coeffs_old = propagate_kinetic_analytical(
            dt / 2, self.params_old, self.coeffs_old, step_context.ones_D
        )
        if self.params_frozen is not None:
            self.params_frozen = self.params_old[step_context.n_dynamic :]

        coeffs_startguess = jnp.zeros(start_guess.shape[0], dtype=self.coeffs_old.dtype)
        start_guess, _ = propagate_kinetic_analytical(dt / 2, start_guess, coeffs_startguess, step_context.ones_D)

        params_solved_temp, coeffs_new_temp, final_RE, nit = self.find_next_timestep_solution(
            start_guess, maxiter=maxiter
        )

        if self.params_frozen is not None:
            params_solved_full = jnp.concatenate((params_solved_temp, self.params_frozen), axis=0)
            params_solved_full, coeffs_new = propagate_kinetic_analytical(
                dt / 2, params_solved_full, coeffs_new_temp, step_context.ones_D
            )
            params_solved = params_solved_full[: params_solved_temp.shape[0]]
            params_frozen_next = params_solved_full[params_solved_temp.shape[0] :]
        else:
            params_solved, coeffs_new = propagate_kinetic_analytical(
                dt / 2, params_solved_temp, coeffs_new_temp, step_context.ones_D
            )

        self.params_old = params_old_temp
        self.coeffs_old = coeffs_old_temp
        self.params_frozen = params_frozen_next

        return StepResult(
            params_solved=params_solved, coeffs_new=coeffs_new, rothe_error=float(final_RE), n_iterations=nit
        )

    def merge_new_and_frozen(self, params_solved):
        if self.params_frozen is not None:
            return jnp.concatenate((params_solved, self.params_frozen), axis=0)
        return params_solved

    def renormalize_coeffs(self, coeffs_new, params_new_full):
        S, H, H2 = self.SHH2(
            t=self.t, params_ket=params_new_full, params_bra=params_new_full, splitting_type=self.splitting_type
        )
        norm = jnp.conj(coeffs_new) @ S @ coeffs_new
        sqrt_norm = jnp.sqrt(norm)
        coeffs_new = coeffs_new / sqrt_norm
        if self.dt.imag != 0:
            print("Norm before renormalization: ", norm)
            new_energy = jnp.real(jnp.conj(coeffs_new) @ H @ coeffs_new)
            H2_expec = jnp.real(jnp.conj(coeffs_new) @ H2 @ coeffs_new)
            print("Energy after renormalization: ", new_energy)
            variance = H2_expec - new_energy**2
            print("Variance after renormalization: ", variance)
        return coeffs_new, S

    def compute_max_overlap(self, S, params_new_full, step_context):
        S_np = np.asarray(S)
        S_abs = np.abs(S_np)
        n_dynamic = params_new_full.shape[0] - step_context.n_frozen
        max_overlap_dyn_dyn = float(np.max(np.triu(S_abs[:n_dynamic, :n_dynamic], k=1))) if n_dynamic > 1 else 0.0
        max_overlap_dyn_frozen = float(np.max(S_abs[:n_dynamic, n_dynamic:])) if step_context.n_frozen > 0 else 0.0
        return max(max_overlap_dyn_dyn, max_overlap_dyn_frozen)

    def compute_observables(self, params_new_full, coeffs_new):
        dipole_moment = None
        double_autocorrelation = None
        if self._dipole_solver is None:
            return dipole_moment, double_autocorrelation

        try:
            params_bra_double = _flip_bra_phase_params(params_new_full)
            self._dipole_solver.update_parameters(params_new_full, params_bra_double)
            S_double = self._dipole_solver.calculate_S()
            double_autocorrelation = jnp.dot(coeffs_new, S_double @ coeffs_new)

            self._dipole_solver.update_parameters(params_new_full, params_new_full)
            dipole_moment = self._dipole_solver.calculate_dipole_moment(coeffs_new)
        except Exception as e:
            print(f"Warning: Could not compute dipole/autocorrelation moment: {e}")

        return dipole_moment, double_autocorrelation

    def save_step_state(self, rothe_error, params_new_full, coeffs_new, step_metrics, n_dynamic):
        if self.output_config.name is None:
            return
        save_rothe_state(
            self.output_config.name,
            self.splitting_type,
            self.output_config.polynomial_string or "",
            self.epsilon,
            self.regularization_lambda,
            self.dt,
            self.t,
            rothe_error,
            params_new_full,
            coeffs_new,
            n_dynamic=n_dynamic,
            rng_state=self._rng.bit_generator.state,
            dipole_moment=step_metrics.dipole_moment,
            double_autocorrelation=step_metrics.double_autocorrelation,
            compression=self.output_config.compression,
            compression_opts=self.output_config.compression_opts,
            path=self.output_config.out_dir,
        )

    def _sample_from_guess_distribution(self, vals, n_samples, n_grid=4096):
        vals = np.asarray(vals, dtype=float)
        vmin = float(np.min(vals))
        vmax = float(np.max(vals))

        x = np.linspace(vmin, vmax, int(n_grid), dtype=float)
        pdf = np.zeros_like(x)
        for center in vals:
            sigma = max(abs(float(center)), 1e-6)
            z = (x - float(center)) / sigma
            pdf += np.exp(-0.5 * z * z) / (np.sqrt(2.0 * np.pi) * sigma)
        pdf_sum = np.sum(pdf)
        if not np.isfinite(pdf_sum) or pdf_sum <= 0.0:
            return self._rng.choice(vals, size=n_samples, replace=True)
        p = pdf / pdf_sum
        return self._rng.choice(x, p=p, size=n_samples)

    def sample_gaussian_candidates(self, params_reference, n_candidates=1000):
        """Sample candidate Gaussians from per-parameter empirical guess distributions."""
        params_reference = np.asarray(params_reference, dtype=float)
        n_ref, D, n_types = params_reference.shape

        candidates = np.zeros((int(n_candidates), D, 4), dtype=float)
        for d in range(D):
            for k in range(4):
                vals = params_reference[:, d, k]
                candidates[:, d, k] = self._sample_from_guess_distribution(vals, int(n_candidates))

        candidates[:, :, 0] = np.maximum(np.abs(candidates[:, :, 0]), 1e-6)
        return candidates

    def candidate_max_overlap(self, params_new_base, candidate_param, t_mid, splitting_eff):
        """Return max |S_ij| between a suggested candidate and the existing basis."""
        candidate_param = jnp.asarray(candidate_param)
        params_existing = jnp.asarray(params_new_base)
        if self.params_frozen is not None:
            params_existing = jnp.concatenate([params_existing, self.params_frozen], axis=0)
        S_candidate, _, _ = self.SHH2(
            t=t_mid, params_ket=params_existing, params_bra=candidate_param, splitting_type=splitting_eff
        )
        return float(np.max(np.abs(np.asarray(S_candidate))))

    def filter_candidate_gaussians_by_overlap(self, params_new_base, candidates, t_mid, splitting_eff):
        """Discard suggested candidates with max overlap above the hard threshold."""
        candidates = np.asarray(candidates, dtype=float)
        if candidates.ndim != 3 or candidates.shape[0] == 0:
            return candidates, 0, np.empty((0,), dtype=float)

        params_existing = jnp.asarray(params_new_base)
        if self.params_frozen is not None:
            params_existing = jnp.concatenate([params_existing, self.params_frozen], axis=0)

        S_candidates, _, _ = self.SHH2(
            t=t_mid, params_ket=params_existing, params_bra=jnp.asarray(candidates), splitting_type=splitting_eff
        )
        S_abs = np.abs(np.asarray(S_candidates))
        if S_abs.shape[0] == candidates.shape[0]:
            max_overlaps = np.max(S_abs, axis=1)
        elif S_abs.shape[-1] == candidates.shape[0]:
            max_overlaps = np.max(S_abs, axis=0)
        else:
            raise ValueError(f"Unexpected overlap matrix shape {S_abs.shape} for {candidates.shape[0]} candidates.")
        keep_mask = max_overlaps <= _CANDIDATE_OVERLAP_DISCARD_THRESHOLD
        kept_candidates = candidates[keep_mask]
        discarded_count = int(candidates.shape[0] - kept_candidates.shape[0])
        return kept_candidates, discarded_count, max_overlaps

    def _compute_base_rothe_blocks(self, params_new_base, SHH2_oldold):
        t_mid = self.t + self.dt / 2
        splitting_eff = "none" if self.splitting_type in (None, "none") else self.splitting_type

        S_xx, H_xx, H2_xx = SHH2_oldold
        dt_real = float(np.real(self.dt))
        dt_imag = float(np.imag(self.dt))
        dt_abs_sq_4 = float(np.abs(self.dt) ** 2 / 4.0)

        B_dagger_B = np.asarray(S_xx + dt_abs_sq_4 * H2_xx + dt_imag * H_xx)

        S_nn, H_nn, H2_nn = self.SHH2(
            t=t_mid, params_ket=params_new_base, params_bra=params_new_base, splitting_type=splitting_eff
        )
        A_dd = np.asarray(S_nn + dt_abs_sq_4 * H2_nn - dt_imag * H_nn)

        S_xn, H_xn, H2_xn = self.SHH2(
            t=t_mid, params_ket=params_new_base, params_bra=self.params_old, splitting_type=splitting_eff
        )
        rho_xd = np.asarray(S_xn - dt_abs_sq_4 * H2_xn - 1j * dt_real * H_xn)

        if self.params_frozen is None:
            A_base = A_dd
            rho_base = rho_xd
            frozen_blocks = None
            S_base = np.asarray(S_nn)
        else:
            n_dyn = params_new_base.shape[0]
            S_ff = np.asarray(S_xx[n_dyn:, n_dyn:])
            H_ff = np.asarray(H_xx[n_dyn:, n_dyn:])
            H2_ff = np.asarray(H2_xx[n_dyn:, n_dyn:])
            A_ff = S_ff + dt_abs_sq_4 * H2_ff - dt_imag * H_ff

            S_df, H_df, H2_df = self.SHH2(
                t=t_mid, params_ket=self.params_frozen, params_bra=params_new_base, splitting_type=splitting_eff
            )
            A_df = np.asarray(S_df + dt_abs_sq_4 * H2_df - dt_imag * H_df)

            S_xf = np.asarray(S_xx[:, n_dyn:])
            H_xf = np.asarray(H_xx[:, n_dyn:])
            H2_xf = np.asarray(H2_xx[:, n_dyn:])
            rho_xf = S_xf - dt_abs_sq_4 * H2_xf - 1j * dt_real * H_xf

            A_base = np.block([[A_dd, A_df], [np.conj(A_df).T, A_ff]])
            rho_base = np.concatenate([rho_xd, rho_xf], axis=1)
            S_base = np.asarray(S_nn)  # dynamic-dynamic only
            frozen_blocks = {"A_df": A_df, "A_ff": A_ff, "rho_xf": rho_xf, "S_df": np.asarray(S_df), "S_ff": S_ff}

        return {
            "A_base": A_base,
            "rho_base": rho_base,
            "B": B_dagger_B,
            "dt_real": dt_real,
            "dt_imag": dt_imag,
            "dt_abs_sq_4": dt_abs_sq_4,
            "t_mid": t_mid,
            "splitting_eff": splitting_eff,
            "frozen_blocks": frozen_blocks,
            "S_base": S_base,
        }

    def _solve_rothe_from_blocks(self, A_dagger_A, rho_mat, B_dagger_B):
        coeffs_old_np = np.asarray(self.coeffs_old)
        rho_vec = np.conj(rho_mat).T @ coeffs_old_np
        A_reg = A_dagger_A + np.eye(A_dagger_A.shape[0]) * self.regularization_lambda
        coeffs_new = np.linalg.solve(A_reg, rho_vec)
        overlap_term = np.conj(coeffs_old_np) @ B_dagger_B @ coeffs_old_np
        projection_term = np.conj(rho_vec).T @ coeffs_new
        rothe_error_val = float(np.real(overlap_term - projection_term))
        return rothe_error_val, jnp.asarray(coeffs_new)

    def _evaluate_candidate_added_error(self, params_new_base, blocks, candidate_param, overlap_penalty_scale=0.0):
        candidate_param = jnp.asarray(candidate_param)
        if candidate_param.ndim != 3 or candidate_param.shape[0] != 1:
            raise ValueError(f"Expected candidate shape (1, D, 4), got {candidate_param.shape}")
        frozen_blocks = blocks["frozen_blocks"]
        has_frozen = frozen_blocks is not None
        A_aug, rho_aug, S_dyn_aug, S_of_aug = _augment_blocks_with_candidate(
            self.SHH2,
            blocks["t_mid"],
            blocks["dt_real"],
            blocks["dt_imag"],
            blocks["dt_abs_sq_4"],
            blocks["splitting_eff"],
            jnp.asarray(blocks["A_base"]),
            jnp.asarray(blocks["rho_base"]),
            jnp.asarray(params_new_base),
            self.params_old,
            candidate_param,
            has_frozen,
            params_frozen=self.params_frozen,
            A_df=jnp.asarray(frozen_blocks["A_df"]) if has_frozen else None,
            A_ff=jnp.asarray(frozen_blocks["A_ff"]) if has_frozen else None,
            rho_xf=jnp.asarray(frozen_blocks["rho_xf"]) if has_frozen else None,
            S_base=jnp.asarray(blocks["S_base"]),
            S_df=jnp.asarray(frozen_blocks["S_df"]) if has_frozen else None,
        )
        err, coeffs_new = self._solve_rothe_from_blocks(np.asarray(A_aug), np.asarray(rho_aug), blocks["B"])
        if overlap_penalty_scale != 0.0:
            S_dyn_frozen_aug = S_of_aug if has_frozen else None
            if has_frozen:
                S_ff = jnp.asarray(frozen_blocks["S_ff"])
                S_new_full = jnp.block([[S_dyn_aug, S_of_aug], [jnp.conj(S_of_aug).T, S_ff]])
            else:
                S_new_full = S_dyn_aug
            # pen = float(np.asarray(_overlap_penalty(S_dyn_aug, S_dyn_frozen_aug) + _log_det_penalty(S_new_full)))
            err = err  # + overlap_penalty_scale * pen
        return err, coeffs_new

    def optimize_single_gaussian(self, params_new_base, candidate_init, SHH2_oldold, maxiter=100):
        """Optimize only the newly added Gaussian while keeping all existing ones fixed."""
        params_new_base = jnp.asarray(params_new_base)
        candidate_init = jnp.asarray(candidate_init)

        if candidate_init.shape[0] != 1:
            raise ValueError(f"Expected candidate_init shape (1, D, 4), got {candidate_init.shape}")

        x0 = np.asarray(candidate_init).ravel()

        def _pack_params(x_flat):
            cand = jnp.asarray(x_flat).reshape(candidate_init.shape)
            cand = cand.at[:, :, 0].set(jnp.maximum(jnp.abs(cand[:, :, 0]), 1e-6))
            return jnp.concatenate([params_new_base, cand], axis=0)

        def f_and_g_single(x_flat):
            params_aug = _pack_params(x_flat)
            val, g = self.rothe_grad_fn(
                params_aug,
                self.params_old,
                self.coeffs_old,
                self.SHH2,
                self.t,
                self.dt,
                self.params_frozen,
                lambda_=self.regularization_lambda,
                SHH2_oldold=SHH2_oldold,
            )
            grad_last = np.asarray(g[-1]).ravel()
            return float(val), grad_last

        if int(maxiter) == 0:
            x_opt = x0
            nit = -1
        else:
            sol = minimize(
                f_and_g_single,
                x0,
                jac=True,
                method="BFGS",
                options={"gtol": 1e-10, "maxiter": int(maxiter), "disp": False},
            )
            x_opt = sol.x
            nit = int(sol.nit)

        params_aug_opt = _pack_params(x_opt)
        final_err, coeffs_new = self._eval_rothe_nograd(params_aug_opt, return_cnew=True, SHH2_oldold=SHH2_oldold)
        return params_aug_opt, coeffs_new, float(final_err), nit

    def reoptimize_all_gaussians(self, params_with_candidate_init, maxiter=100):
        """Re-optimize all dynamic Gaussians after adding a new one."""
        params_solved, coeffs_new, final_RE, nit = self.find_next_timestep_solution(
            params_init=np.asarray(params_with_candidate_init), maxiter=maxiter, given_gtol=1e-14
        )
        return params_solved, coeffs_new, float(final_RE), int(nit)

    def adding_strategy_as_a_whole(
        self, params_new_base, coeffs_before, rothe_error_before, nit_before, threshold, maxiter
    ):
        """Sample candidates, add best Gaussian, optimize it, then re-optimize all Gaussians."""
        print(f"A Gaussian is added because the Rothe error is bigger than {threshold:.3e}.")
        print(f"Squared Rothe error before adding a Gaussian: {rothe_error_before:.3e}")

        SHH2_oldold = compute_SHH2_oldold(self.SHH2, self.params_old, self.t + self.dt / 2, self.splitting_type)
        base_blocks = self._compute_base_rothe_blocks(params_new_base, SHH2_oldold)

        params_dynamic_only = np.asarray(params_new_base)
        candidates = self.sample_gaussian_candidates(params_dynamic_only, n_candidates=1000)
        candidates, discarded_overlap_count, _ = self.filter_candidate_gaussians_by_overlap(
            params_new_base, candidates, t_mid=base_blocks["t_mid"], splitting_eff=base_blocks["splitting_eff"]
        )
        if discarded_overlap_count > 0:
            print(
                f"Discarded {discarded_overlap_count} candidate Gaussian(s) with max overlap "
                f"> {_CANDIDATE_OVERLAP_DISCARD_THRESHOLD:.2f}."
            )
        print(f"Evaluating {candidates.shape[0]} candidate Gaussian(s) after overlap filtering.")

        if candidates.shape[0] == 0:
            print("No admissible Gaussian candidates remain after overlap filtering; skipping Gaussian addition.")
            return StepResult(
                params_solved=jnp.asarray(params_new_base),
                coeffs_new=jnp.asarray(coeffs_before),
                rothe_error=float(rothe_error_before),
                n_iterations=int(nit_before),
            )

        best_idx = None
        best_err = np.inf
        best_candidate = None
        for idx in range(candidates.shape[0]):
            candidate = jnp.asarray(candidates[idx : idx + 1])
            err, _ = self._evaluate_candidate_added_error(
                params_new_base, base_blocks, candidate, overlap_penalty_scale=float(0)
            )
            if err < best_err:
                best_err = float(err)
                best_idx = idx
                best_candidate = candidate

        if best_idx is None or best_candidate is None:
            raise RuntimeError(
                "Could not find a finite Rothe-error Gaussian candidate after overlap-based candidate discarding."
            )

        # Compute pure Rothe error (without overlap penalty) for the best candidate
        best_err_pure, _ = self._evaluate_candidate_added_error(
            params_new_base, base_blocks, best_candidate, overlap_penalty_scale=0.0
        )
        print(f"Squared Rothe error after adding a Gaussian: {best_err_pure:.3e} (penalized: {best_err:.3e})")

        params_after_single_opt, coeffs_after_single_opt, err_single_opt, _ = self.optimize_single_gaussian(
            params_new_base, best_candidate, SHH2_oldold=SHH2_oldold, maxiter=max(20, int(maxiter // 2))
        )
        print(f"Squared Rothe error after optimizing one Gaussian: {err_single_opt:.3e}")

        params_after_full_reopt, coeffs_after_full_reopt, err_full_reopt, nit_full = self.reoptimize_all_gaussians(
            params_after_single_opt, maxiter=max(int(maxiter), 1000)
        )
        print(f"Squared Rothe error after full reoptimization: {err_full_reopt:.3e}")

        return StepResult(
            params_solved=params_after_full_reopt,
            coeffs_new=coeffs_after_full_reopt,
            rothe_error=float(err_full_reopt),
            n_iterations=int(nit_full),
        )

    def advance_solver_state(self, params_new_full, coeffs_new, start_time):
        n_frozen = self.params_frozen.shape[0] if self.params_frozen is not None else 0
        n_dynamic = self.params_old.shape[0] - n_frozen
        self.params_oldold = self.params_old[:n_dynamic]  # dynamic part only
        self.params_old = params_new_full
        self.coeffs_old = coeffs_new
        time_taken = time.time() - start_time
        self.total_time += time_taken
        print(f"Time taken: {time_taken:.2f}")

    def propagate(
        self, num_iterations, num_time_steps_total, params_init=None, coeffs_init=None, maxiter=100, renormalize=True
    ):
        self.init_step_context(params_init=params_init, coeffs_init=coeffs_init)
        self.threshold = self.epsilon / float(num_time_steps_total)

        for _ in range(num_iterations):
            added_gaussian_this_step = False
            step_context = self.init_step_context()
            self.t += self.dt
            start = time.time()
            start_guess = self.build_startguess(step_context)
            if self.splitting_type in (None, "none"):
                step_result = self.solve_step_none(start_guess, maxiter=maxiter)
            elif self.splitting_type == "kinetic":
                step_result = self.solve_step_kinetic(start_guess, step_context, maxiter=maxiter)
            else:
                raise ValueError(f"Unknown splitting_type {self.splitting_type!r}; expected 'none' or 'kinetic'.")
            if step_result.rothe_error < 0:
                print(
                    f"Warning: Negative squared Rothe error {step_result.rothe_error:.3e} encountered. "
                    "This may indicate numerical issues. Setting error to 0 for threshold comparison."
                )
                sqrt_re = 0.0
            else:
                sqrt_re = float(step_result.rothe_error) ** 0.5
            if sqrt_re > self.threshold:
                print(
                    f"Squared Rothe error too big ({step_result.rothe_error:.3e}), "
                    f"rerunning optimization with maxiter=1000."
                )
                if self.splitting_type in (None, "none"):
                    step_result = self.solve_step_none(step_result.params_solved, maxiter=1000)
                elif self.splitting_type == "kinetic":
                    step_result = self.solve_step_kinetic(step_result.params_solved, step_context, maxiter=1000)
                print(f"Squared Rothe error after rerun: {step_result.rothe_error:.3e}")

                if float(step_result.rothe_error) ** 0.5 > self.threshold:
                    if self._gaussian_addition_cooldown_remaining > 0:
                        print(
                            "Skipping Gaussian addition due to cooldown "
                            f"({self._gaussian_addition_cooldown_remaining} step(s) remaining)."
                        )
                    else:
                        step_result = self.adding_strategy_as_a_whole(
                            params_new_base=step_result.params_solved,
                            coeffs_before=step_result.coeffs_new,
                            rothe_error_before=float(step_result.rothe_error),
                            nit_before=int(step_result.n_iterations),
                            threshold=self.threshold,
                            maxiter=maxiter,
                        )
                        added_gaussian_this_step = True
                        self._gaussian_addition_cooldown_remaining = self.gaussian_addition_cooldown_steps

            params_new_full = self.merge_new_and_frozen(step_result.params_solved)
            if np.isrealobj(self.t) and not np.iscomplexobj(self.t):
                time_str = f"{float(np.real(self.t)):.2f}"
            else:
                time_str = str(self.t)

            coeffs_new = step_result.coeffs_new
            max_overlap = float("nan")
            if renormalize:
                coeffs_new, S = self.renormalize_coeffs(coeffs_new, params_new_full)
                max_overlap = self.compute_max_overlap(S, params_new_full, step_context)

            self.cumulative_rothe_error += abs(float(step_result.rothe_error)) ** 0.5
            print(
                f"Time {time_str}:Squared Rothe error: {step_result.rothe_error:.3e}, "
                f"Number of iterations: {step_result.n_iterations}, "
                f"Max overlap: {max_overlap:.4f}, "
                f"Cumulative Rothe error: {self.cumulative_rothe_error:.2e}"
            )

            dipole_moment, double_autocorrelation = self.compute_observables(params_new_full, coeffs_new)
            step_metrics = StepMetrics(
                max_overlap=max_overlap, dipole_moment=dipole_moment, double_autocorrelation=double_autocorrelation
            )
            self.save_step_state(
                step_result.rothe_error,
                params_new_full,
                coeffs_new,
                step_metrics,
                n_dynamic=int(step_result.params_solved.shape[0]),
            )
            self.previous_step_rothe_error = float(step_result.rothe_error)
            self.advance_solver_state(params_new_full, coeffs_new, start)
            if not added_gaussian_this_step and self._gaussian_addition_cooldown_remaining > 0:
                self._gaussian_addition_cooldown_remaining -= 1
        print(f"Total time: {self.total_time}")

    def resume_from_file(self, t_start=None):
        """
        Load state from HDF5 and set solver to resume at requested time.

        If t_start is None -> resume from last saved step.
        If t_start exists -> resume at that step and delete later steps.
        If t_start doesn't exist -> resume from latest step < t_start.
        """
        if self.output_config.name is None:
            raise ValueError("Cannot resume: 'name' is None (no filename to load).")
        filename = os.path.join(self.output_config.out_dir, f"{self.output_config.name}__{self.splitting_type}.h5")
        saved = load_rothe_file(self.output_config.name, self.splitting_type, path=self.output_config.out_dir)
        saved_lambda = float(saved["attrs"]["regularization_lambda"])
        if not np.isclose(saved_lambda, self.regularization_lambda, rtol=0.0, atol=0.0):
            raise ValueError(
                "Requested regularization_lambda does not match the saved file: "
                f"{self.regularization_lambda} != {saved_lambda}"
            )
        sel = _select_and_trim_to_time(filename, t_start)
        params_saved = np.asarray(sel["params"])
        coeffs_saved = np.asarray(sel["coeffs"])

        if params_saved.ndim != self.params_old.ndim or params_saved.shape[1:] != tuple(self.params_old.shape[1:]):
            raise ValueError(
                f"Saved params tail shape {params_saved.shape[1:]} != current {tuple(self.params_old.shape[1:])}"
            )
        if coeffs_saved.shape[0] != params_saved.shape[0]:
            raise ValueError(f"Saved coeffs size {coeffs_saved.shape[0]} != saved params count {params_saved.shape[0]}")

        # Use the full saved size — may be larger than initial if a Gaussian was added mid-run.
        self.params_old = jnp.asarray(params_saved)
        self.coeffs_old = jnp.asarray(coeffs_saved)
        self.t = np.asarray(sel["t"]).item()

        # Reconstruct params_oldold as dynamic-only from the previous step.
        if sel["params_prev"] is not None and sel["n_dynamic_prev"] is not None:
            n_dyn_prev = sel["n_dynamic_prev"]
            self.params_oldold = jnp.asarray(sel["params_prev"][:n_dyn_prev])
        else:
            self.params_oldold = None

        self.previous_step_rothe_error = float(saved["rothe_error"][sel["idx"]])

        rng_state = sel.get("rng_state")
        if isinstance(rng_state, dict):
            self._rng.bit_generator.state = rng_state

        self.cumulative_rothe_error = float(np.sum(np.sqrt(np.abs(saved["rothe_error"]))))
        return sel
