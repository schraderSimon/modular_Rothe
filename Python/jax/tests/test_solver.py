"""Tests for the RotheSolver class (rothe/solver.py).

All tests use the 1D quantum harmonic oscillator (ω=1) whose ground state is
exactly representable by a single Gaussian.  This makes it possible to assess:

  1. convergence of the Rothe error after BFGS optimisation,
  2. norm conservation during real-time propagation,
  3. consistency of the Rothe error with the expected dt scaling.
"""

import numpy as np
import pytest

import jax.numpy as jnp
from rothe.solver import RotheSolver, setUpRotheErrorAndGradient_jit
from rothe.systems.common import set_up_SHH2

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_POLY_STRING_1D_HO = "dimension 1\npolynomial\nx0x0: 0.5\n"


def _ground_state_params_1d():
    """Single Gaussian that is the exact 1D HO ground state (a = 1/√2)."""
    params = jnp.zeros((1, 1, 4))
    params = params.at[0, 0, 0].set(1.0 / jnp.sqrt(2.0))
    return params


def _make_solver(dt):
    """Return a fresh RotheSolver + SHH2 for the 1D HO ground state."""
    params = _ground_state_params_1d()
    coeffs = jnp.array([1.0 + 0.0j])
    SHH2 = set_up_SHH2(potential_string=_POLY_STRING_1D_HO, D=1)
    rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit("none")
    solver = RotheSolver(
        SHH2=SHH2,
        dt=dt,
        t=0.0,
        epsilon=1.0,
        regularization_lambda=1e-12,
        params_old=params,
        coeffs_old=coeffs,
        rothe_grad_fn=rothe_vg_jit,
        rothe_nograd=rothe_error,
    )
    return solver, SHH2


# ---------------------------------------------------------------------------
# Test 1 – Rothe error is negligible after BFGS optimisation
# ---------------------------------------------------------------------------


def test_rothe_error_small_after_optimization():
    """After optimisation the Rothe error for the HO ground state must be tiny.

    For an exact eigenstate the Rothe functional has a global minimum of zero
    (achieved by keeping the same Gaussian and rotating the phase of c).
    """
    solver, _ = _make_solver(dt=0.05)
    _, _, final_RE, _ = solver.find_next_timestep_solution(maxiter=50)
    assert final_RE < 1e-6, f"Optimised Rothe error too large: {final_RE:.3e}"


# ---------------------------------------------------------------------------
# Test 2 – Norm is conserved after a single propagation step
# ---------------------------------------------------------------------------


def test_norm_conserved_after_single_step():
    """Without explicit renormalisation the overlap ⟨ψ|ψ⟩ should stay ≈ 1.

    Real-time propagation by a unitary operator preserves norm exactly; the
    discrete Rothe approximation should recover this up to a small tolerance.
    """
    solver, SHH2 = _make_solver(dt=0.05)
    params_new, coeffs_new, _, _ = solver.find_next_timestep_solution(maxiter=50)
    S, _, _ = SHH2(t=solver.dt, params_ket=params_new, params_bra=params_new)
    norm = float(jnp.real(jnp.vdot(coeffs_new, S @ coeffs_new)))
    assert abs(norm - 1.0) < 1e-4, f"Norm after one step: {norm:.8f}"


# ---------------------------------------------------------------------------
# Test 3 – Unoptimised Rothe error decreases with finer dt (consistency)
# ---------------------------------------------------------------------------


def test_rothe_error_decreases_with_finer_dt():
    """The identity-guess Rothe error must be smaller for a smaller time step.

    A displaced coherent state (mu=1.0) has non-trivial real-time dynamics, so
    the identity-guess Rothe error is appreciable and shrinks clearly as dt
    decreases — confirming consistency of the time-stepping scheme.

    The exact eigenstate cannot be used here because its Rothe error is already
    at machine-precision for any dt, making the comparison meaningless.
    """
    # Displaced coherent state: same width as ground state but centered at mu=1
    params_displaced = jnp.zeros((1, 1, 4))
    params_displaced = params_displaced.at[0, 0, 0].set(1.0 / jnp.sqrt(2.0))
    params_displaced = params_displaced.at[0, 0, 2].set(1.0)  # mu = 1
    coeffs = jnp.array([1.0 + 0.0j])

    def _make_displaced_solver(dt):
        SHH2 = set_up_SHH2(potential_string=_POLY_STRING_1D_HO, D=1)
        rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit("none")
        return RotheSolver(
            SHH2=SHH2,
            dt=dt,
            t=0.0,
            epsilon=1.0,
            regularization_lambda=1e-12,
            params_old=params_displaced,
            coeffs_old=coeffs,
            rothe_grad_fn=rothe_vg_jit,
            rothe_nograd=rothe_error,
        )

    _, _, err_coarse, _ = _make_displaced_solver(dt=0.2).find_next_timestep_solution(maxiter=0)
    _, _, err_fine, _ = _make_displaced_solver(dt=0.1).find_next_timestep_solution(maxiter=0)
    assert err_fine < err_coarse, f"Expected err(dt=0.1)={err_fine:.3e} < err(dt=0.2)={err_coarse:.3e}"


# ---------------------------------------------------------------------------
# Test 4 – Norm is conserved over multiple propagation steps
# ---------------------------------------------------------------------------


def test_norm_conserved_over_multiple_steps():
    """Propagate for several steps without renormalisation; norm must stay ≈ 1.

    Uses RotheSolver.propagate() with renormalize=False so that norm stays
    are not enforced by the solver itself.
    """
    dt = 0.05
    n_steps = 4

    solver, SHH2 = _make_solver(dt=dt)
    solver.propagate(num_iterations=n_steps, num_time_steps_total=n_steps, maxiter=50, renormalize=False)

    # solver.params_old and solver.coeffs_old now hold the state after n_steps
    S, _, _ = SHH2(t=float(n_steps * dt), params_ket=solver.params_old, params_bra=solver.params_old)
    norm = float(jnp.real(jnp.vdot(solver.coeffs_old, S @ solver.coeffs_old)))
    assert abs(norm - 1.0) < 2e-4, f"Norm after {n_steps} steps: {norm:.8f} (expected ≈ 1.0)"


# ---------------------------------------------------------------------------
# Test 5 – Energy is conserved over multiple propagation steps
# ---------------------------------------------------------------------------


def test_energy_conserved_over_multiple_steps():
    """The expectation value ⟨H⟩ should remain ≈ E₀ = 0.5 throughout propagation."""
    dt = 0.05
    n_steps = 4
    E0_expected = 0.5  # ground state of 1D HO with ω=1

    solver, SHH2 = _make_solver(dt=dt)
    solver.propagate(num_iterations=n_steps, num_time_steps_total=n_steps, maxiter=50, renormalize=True)

    S, H, _ = SHH2(t=float(n_steps * dt), params_ket=solver.params_old, params_bra=solver.params_old)
    # Renormalise coefficient vector before computing expectation value
    norm = jnp.real(jnp.vdot(solver.coeffs_old, S @ solver.coeffs_old))
    c = solver.coeffs_old / jnp.sqrt(norm)
    energy = float(jnp.real(jnp.vdot(c, H @ c)))
    assert abs(energy - E0_expected) < 1e-3, f"Energy after {n_steps} steps: {energy:.6f} (expected {E0_expected})"


def test_candidate_overlap_detects_duplicate_gaussian():
    """Suggested candidates identical to an existing Gaussian should be discarded."""
    solver, _ = _make_solver(dt=0.05)
    params_base = _ground_state_params_1d()
    candidate_duplicate = params_base[:1]
    max_overlap = solver.candidate_max_overlap(
        params_base, candidate_duplicate, t_mid=solver.t + solver.dt / 2, splitting_eff="none"
    )
    assert max_overlap > 0.99, f"Expected duplicate candidate overlap > 0.99, got {max_overlap:.6f}"


def test_overlap_filter_discards_duplicate_candidates():
    """Batch overlap screening should remove duplicate suggested Gaussians."""
    solver, _ = _make_solver(dt=0.05)
    params_base = _ground_state_params_1d()
    candidates = jnp.concatenate([params_base[:1], params_base[:1]], axis=0)
    kept_candidates, discarded_count, max_overlaps = solver.filter_candidate_gaussians_by_overlap(
        params_base, candidates, t_mid=solver.t + solver.dt / 2, splitting_eff="none"
    )
    assert discarded_count == 2, f"Expected 2 discarded candidates, got {discarded_count}"
    assert kept_candidates.shape[0] == 0, f"Expected 0 kept candidates, got {kept_candidates.shape[0]}"
    assert np.all(max_overlaps > 0.99), f"Expected all max overlaps > 0.99, got {max_overlaps}"
