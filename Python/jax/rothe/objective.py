"""
Rothe objective assembly and JAX-compiled objective/gradient builders.

This module contains the low-level block-matrix assembly used by the Rothe
objective and the factory for creating JIT-compiled value/gradient callables.
"""

from dataclasses import dataclass
from functools import partial

import jax
import jax.numpy as jnp
from jax.scipy.linalg import solve
from rothe.optimization import RotheObjectiveContext


@dataclass(frozen=True)
class BasisParams:
    """Compact container for basis-parameter blocks used in objective assembly."""

    new: jnp.ndarray
    old: jnp.ndarray
    frozen: jnp.ndarray | None = None


def _overlap_penalty(S_dyn_dyn, S_dyn_frozen=None, center=0.97**2, sharpness=100.0, multiplier=1.0):
    """Softplus penalty for near-linearly-dependent Gaussian pairs.

    Uses a shifted softplus (smooth ReLU) in |S_ij|²:
        pen(x) = softplus(sharpness * (x - center)) / sharpness
    This is effectively zero for |S_ij|² << center and grows linearly
    for |S_ij|² >> center, with a smooth transition of width ~1/sharpness.

    At center=0.97, sharpness=100 the transition width is ~0.01 in |S_ij|²,
    so the penalty is negligible below ~0.96 and clearly active above ~0.98.

    Dynamic-dynamic pairs are counted once (upper triangle).
    Dynamic-frozen pairs (when S_dyn_frozen is not None) are all counted.
    """

    def _softplus_pen(mat):
        abs2 = jnp.abs(mat) ** 2
        return jax.nn.softplus(sharpness * (abs2 - center)) / sharpness

    pen = jnp.sum(jnp.triu(_softplus_pen(S_dyn_dyn), k=1))
    if S_dyn_frozen is not None:
        pen = pen + jnp.sum(_softplus_pen(S_dyn_frozen))
    return pen


def _log_det_penalty(S_new_full, multiplier=100):
    """Penalty -log(det(S)) / N on the full new-time basis overlap matrix.

    Discourages near-linear dependence across all new basis functions (dynamic
    and frozen).  As the basis approaches linear dependence, det(S) -> 0 and
    the penalty diverges to +inf, pushing parameters apart.
    """
    S_herm = 0.5 * (S_new_full + jnp.conj(S_new_full).T)
    eigvals = jnp.real(jnp.linalg.eigvalsh(S_herm))
    eigvals_clipped = jnp.clip(eigvals, 1e-14, jnp.inf)
    return -jnp.mean(jnp.log(eigvals_clipped)) / multiplier


@partial(jax.jit, static_argnames=("SHH2", "splitting_type"))
def compute_SHH2_oldold(SHH2, params_old, t_mid, splitting_type):
    """Precompute S, H, H² at midpoint time for old-old block (incl. frozen terms)."""
    return SHH2(t=t_mid, params_ket=params_old, params_bra=params_old, splitting_type=splitting_type)


def build_objective_context(SHH2, params_old, coeffs_old, t, dt, splitting_type, params_frozen):
    """Build objective context and precompute SHH2_oldold for the current solver state."""
    t_mid = t + dt / 2
    SHH2_oldold = compute_SHH2_oldold(SHH2, params_old, t_mid, splitting_type)
    return RotheObjectiveContext(
        params_old=params_old,
        coeffs_old=coeffs_old,
        SHH2=SHH2,
        t=t,
        dt=dt,
        params_frozen=params_frozen,
        SHH2_oldold=SHH2_oldold,
    )


def _get_oldold_blocks(SHH2, params_old, t_mid, splitting_eff, SHH2_oldold):
    if SHH2_oldold is None:
        SHH2_oldold = compute_SHH2_oldold(SHH2, params_old, t_mid, splitting_eff)
    return SHH2_oldold


def _compute_frozen_coupling_blocks(SHH2, t_mid, splitting_eff, basis):
    """Compute SHH2 blocks involving optimizable and frozen basis functions."""
    S_oo, H_oo, H2_oo = SHH2(t=t_mid, params_ket=basis.new, params_bra=basis.new, splitting_type=splitting_eff)
    S_of, H_of, H2_of = SHH2(t=t_mid, params_ket=basis.frozen, params_bra=basis.new, splitting_type=splitting_eff)
    S_xo, H_xo, H2_xo = SHH2(t=t_mid, params_ket=basis.new, params_bra=basis.old, splitting_type=splitting_eff)
    return (S_oo, H_oo, H2_oo), (S_of, H_of, H2_of), (S_xo, H_xo, H2_xo)


def _assemble_AdA_and_rho(SHH2, splitting_eff, t_mid, dt, basis, SHH2_oldold):
    """Assemble linear-system blocks for one Rothe objective evaluation.

    Notation:
    - ``n``: optimizable new basis
    - ``o``: optimizable old basis
    - ``f``: frozen basis
    - ``x``: full old-state basis (``o`` + ``f`` when frozen is present;
      otherwise ``x == o``)

    Shapes:
    - ``basis.new``: ``(n_new, D, 4)``
    - ``basis.old``: ``(n_old, D, 4)``
    - ``SHH2_oldold`` entries: each ``(n_old, n_old)``

    Returns:
    - ``A_dagger_A``: ``(n_eff, n_eff)``
    - ``rho_mat``: ``(n_old, n_eff)``
    - ``B_dagger_B``: ``(n_old, n_old)``

    where ``n_eff = n_new`` if no frozen basis is used, otherwise
    ``n_eff = n_new + n_frozen``.
    """
    dt_real = jnp.real(dt)
    dt_imag = jnp.imag(dt)
    dt_abs_sq_4 = jnp.abs(dt) ** 2 / 4
    S_xx, H_xx, H2_xx = SHH2_oldold
    B_dagger_B = S_xx + dt_abs_sq_4 * H2_xx + dt_imag * H_xx

    if basis.frozen is not None:
        n_opt = basis.old.shape[0] - basis.frozen.shape[0]
        S_ff = S_xx[n_opt:, n_opt:]
        H_ff = H_xx[n_opt:, n_opt:]
        H2_ff = H2_xx[n_opt:, n_opt:]
        AdA_ff = S_ff + dt_abs_sq_4 * H2_ff - dt_imag * H_ff

        S_xf, H_xf, H2_xf = S_xx[:, n_opt:], H_xx[:, n_opt:], H2_xx[:, n_opt:]
        rho_xf = S_xf - dt_abs_sq_4 * H2_xf - 1j * dt_real * H_xf
        oo, of, xo = _compute_frozen_coupling_blocks(SHH2, t_mid, splitting_eff, basis)
        (S_oo, H_oo, H2_oo), (S_of, H_of, H2_of), (S_xo, H_xo, H2_xo) = oo, of, xo

        AdA_oo = S_oo + dt_abs_sq_4 * H2_oo - dt_imag * H_oo
        AdA_of = S_of + dt_abs_sq_4 * H2_of - dt_imag * H_of
        rho_xo = S_xo - dt_abs_sq_4 * H2_xo - 1j * dt_real * H_xo

        A_dagger_A = jnp.block([[AdA_oo, AdA_of], [jnp.conj(AdA_of).T, AdA_ff]])
        rho_mat = jnp.concatenate([rho_xo, rho_xf], axis=1)
        S_new_full = jnp.block([[S_oo, S_of], [jnp.conj(S_of).T, S_ff]])
        penalty_mats = (S_oo, S_of, S_new_full)
    else:
        S_nn, H_nn, H2_nn = SHH2(t=t_mid, params_ket=basis.new, params_bra=basis.new, splitting_type=splitting_eff)
        S_xn, H_xn, H2_xn = SHH2(t=t_mid, params_ket=basis.new, params_bra=basis.old, splitting_type=splitting_eff)
        A_dagger_A = S_nn + dt_abs_sq_4 * H2_nn - dt_imag * H_nn
        rho_mat = S_xn - dt_abs_sq_4 * H2_xn - 1j * dt_real * H_xn
        penalty_mats = (S_nn, None, S_nn)

    return A_dagger_A, rho_mat, B_dagger_B, penalty_mats


def _solve_coeffs_new(A_dagger_A, rho_mat, coeffs_old, lambda_):
    rho_vec = jnp.conj(rho_mat).T @ coeffs_old
    I = jnp.eye(A_dagger_A.shape[0])
    S_reg = A_dagger_A + I * lambda_
    coeffs_new = solve(S_reg, rho_vec)
    return coeffs_new, rho_vec


def setUpRotheErrorAndGradient_jit(splitting_type):
    """
    Build the Rothe error function and its JIT-compiled value-and-gradient variant.

    Returns:
        rothe_error: callable
            ``(params_new, params_old, coeffs_old, SHH2, t, dt,
               params_frozen=None, return_cnew=False, lambda_=1e-8) -> error``

            ``params_frozen`` (optional): nonlinear basis-function parameters that are
            included in the new-time wave function but held fixed during optimisation.
            The full new-time basis is ``concat(params_new, params_frozen)``.
            The gradient is always w.r.t. ``params_new`` only (argnums=0).

        rothe_vg_jit: JIT-compiled value_and_grad of rothe_error w.r.t. params_new
    """

    def rothe_error(
        params_new,
        params_old,
        coeffs_old,
        SHH2,
        t,
        dt,
        params_frozen=None,
        return_cnew=False,
        lambda_=1e-10,
        SHH2_oldold=None,
        overlap_penalty_scale=0.0,
    ):
        splitting_eff = "none" if splitting_type in (None, "none") else splitting_type
        t_mid = t + dt / 2
        basis = BasisParams(new=params_new, old=params_old, frozen=params_frozen)
        SHH2_oldold = _get_oldold_blocks(SHH2, basis.old, t_mid, splitting_eff, SHH2_oldold)
        A_dagger_A, rho_mat, B_dagger_B, penalty_mats = _assemble_AdA_and_rho(
            SHH2, splitting_eff, t_mid, dt, basis, SHH2_oldold
        )
        coeffs_new, rho_vec = _solve_coeffs_new(A_dagger_A, rho_mat, coeffs_old, lambda_)

        overlap_term = jnp.conj(coeffs_old) @ B_dagger_B @ coeffs_old
        projection_term = jnp.conj(rho_vec).T @ coeffs_new
        rothe_error_val = jnp.real(overlap_term - projection_term)

        if not return_cnew:
            S_dyn_dyn, S_dyn_frozen, S_new_full = penalty_mats
            multiplier = 1
            penalty = (
                multiplier * jnp.abs(overlap_penalty_scale) * _overlap_penalty(S_dyn_dyn, S_dyn_frozen, center=0.99**2)
            )
            rothe_error_val = rothe_error_val + penalty

        if return_cnew:
            return rothe_error_val, coeffs_new
        else:
            return rothe_error_val

    rothe_vg_jit = jax.jit(jax.value_and_grad(rothe_error, argnums=0), static_argnames=("SHH2", "return_cnew"))
    rothe_error_jit = jax.jit(rothe_error, static_argnames=("SHH2", "return_cnew"))
    return rothe_error_jit, rothe_vg_jit
