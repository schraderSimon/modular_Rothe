"""
Rothe objective assembly and JAX-compiled objective/gradient builders.

This module contains the low-level block-matrix assembly used by the Rothe
objective and the factory for creating JIT-compiled value/gradient callables.
"""

from dataclasses import dataclass
from functools import partial

import jax
import jax.numpy as jnp
from rothe.block_assembly import make_A_block, make_B_block, make_rho_block, shh2_triplets, split_oldold_frozen_triplets
from rothe.linear_solve import solve_and_score_jax
from rothe.optimization import RotheObjectiveContext


@dataclass(frozen=True)
class BasisParams:
    """Compact container for basis-parameter blocks used in objective assembly."""

    new: jnp.ndarray
    old: jnp.ndarray
    frozen: jnp.ndarray | None = None


def _overlap_penalty(S_dyn_dyn, S_dyn_frozen=None, center=0.99**2, sharpness=100.0, multiplier=1.0):
    """Softplus penalty for near-linearly-dependent Gaussian pairs."""
    norm = jax.nn.softplus(sharpness)

    def _softplus_pen(mat):
        abs2 = jnp.abs(mat) ** 2
        z = (abs2 - center) / (1.0 - center)
        return jax.nn.softplus(sharpness * z) / norm

    pen = jnp.sum(jnp.triu(_softplus_pen(S_dyn_dyn), k=1))

    if S_dyn_frozen is not None:
        pen = pen + jnp.sum(_softplus_pen(S_dyn_frozen))

    return multiplier * pen


@partial(jax.jit, static_argnames=("SHH2", "splitting_type"))
def compute_SHH2_oldold(SHH2, params_old, t_mid, splitting_type):
    """Precompute S, H, H² at midpoint time for old-old block (incl. frozen terms)."""
    return SHH2(t=t_mid, params_ket=params_old, params_bra=params_old, splitting_type=splitting_type)


def build_objective_context(SHH2, params_old, coeffs_old, t, dt, splitting_type, params_frozen, regularization_lambda):
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
        regularization_lambda=regularization_lambda,
    )


def _assemble_AdA_and_rho(SHH2, splitting_eff, t_mid, dt, basis, SHH2_oldold):
    """Assemble linear-system blocks for one Rothe objective evaluation.

    Notation:
    - ``n``: optimizable new basis
    - ``o``: dynamic old basis (optimizable at the previous step)
    - ``f``: frozen basis
        - ``x``: full old-state basis (``o`` + ``f`` when frozen is present,
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
    B_dagger_B = make_B_block(S_xx, H_xx, H2_xx, dt_abs_sq_4, dt_imag)

    if basis.frozen is not None:
        n_dyn = basis.old.shape[0] - basis.frozen.shape[0]
        (S_ff, H_ff, H2_ff), (S_xf, H_xf, H2_xf) = split_oldold_frozen_triplets(S_xx, H_xx, H2_xx, n_dyn)

        AdA_ff = make_A_block(S_ff, H_ff, H2_ff, dt_abs_sq_4, dt_imag)
        rho_xf = make_rho_block(S_xf, H_xf, H2_xf, dt_abs_sq_4, dt_real)

        (S_nn, H_nn, H2_nn), (S_nf, H_nf, H2_nf), (S_xn, H_xn, H2_xn) = shh2_triplets(
            SHH2, t_mid, splitting_eff, ((basis.new, basis.new), (basis.frozen, basis.new), (basis.new, basis.old))
        )

        AdA_nn = make_A_block(S_nn, H_nn, H2_nn, dt_abs_sq_4, dt_imag)
        AdA_nf = make_A_block(S_nf, H_nf, H2_nf, dt_abs_sq_4, dt_imag)
        rho_xn = make_rho_block(S_xn, H_xn, H2_xn, dt_abs_sq_4, dt_real)

        A_dagger_A = jnp.block([[AdA_nn, AdA_nf], [jnp.conj(AdA_nf).T, AdA_ff]])
        rho_mat = jnp.concatenate([rho_xn, rho_xf], axis=1)
        S_new_full = jnp.block([[S_nn, S_nf], [jnp.conj(S_nf).T, S_ff]])
        H_new_full = jnp.block([[H_nn, H_nf], [jnp.conj(H_nf).T, H_ff]])
        H2_new_full = jnp.block([[H2_nn, H2_nf], [jnp.conj(H2_nf).T, H2_ff]])
        penalty_mats = (S_nn, S_nf, S_new_full)
        prediction_mats = (S_new_full, H_new_full, H2_new_full)
    else:
        (S_nn, H_nn, H2_nn), (S_xn, H_xn, H2_xn) = shh2_triplets(
            SHH2, t_mid, splitting_eff, ((basis.new, basis.new), (basis.new, basis.old))
        )
        A_dagger_A = make_A_block(S_nn, H_nn, H2_nn, dt_abs_sq_4, dt_imag)
        rho_mat = make_rho_block(S_xn, H_xn, H2_xn, dt_abs_sq_4, dt_real)
        penalty_mats = (S_nn, None, S_nn)
        prediction_mats = (S_nn, H_nn, H2_nn)

    return A_dagger_A, rho_mat, B_dagger_B, penalty_mats, prediction_mats


def setUpRotheErrorAndGradient_jit(splitting_type):
    """
    Build the Rothe error function and its JIT-compiled value-and-gradient variant.

    Returns:
        rothe_error: callable
            ``(params_new, params_old, coeffs_old, SHH2, t, dt,
               params_frozen=None, return_cnew=False, lambda_) -> error``

            ``lambda_`` must be passed explicitly; no internal default is used.

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
        lambda_=None,
        SHH2_oldold=None,
        overlap_penalty_scale=0.0,
        predictive_penalty_scale=0.0,
    ):
        if lambda_ is None:
            raise ValueError("lambda_ must be provided explicitly (use the user-provided regularization_lambda).")
        splitting_eff = "none" if splitting_type in (None, "none") else splitting_type
        t_mid = t + dt / 2
        basis = BasisParams(new=params_new, old=params_old, frozen=params_frozen)
        if SHH2_oldold is None:
            SHH2_oldold = compute_SHH2_oldold(SHH2, basis.old, t_mid, splitting_eff)
        A_dagger_A, rho_mat, B_dagger_B, penalty_mats, prediction_mats = _assemble_AdA_and_rho(
            SHH2, splitting_eff, t_mid, dt, basis, SHH2_oldold
        )
        rothe_error_val, coeffs_new = solve_and_score_jax(A_dagger_A, rho_mat, B_dagger_B, coeffs_old, lambda_)

        if not return_cnew:
            S_dyn_dyn, S_dyn_frozen, S_new_full = penalty_mats
            S_pred_full, H_pred_full, H2_pred_full = prediction_mats
            multiplier = 0.01
            penalty = (
                multiplier * jnp.abs(overlap_penalty_scale) * _overlap_penalty(S_dyn_dyn, S_dyn_frozen, center=0.990**2)
            )
            rothe_error_val = rothe_error_val + penalty

        if return_cnew:
            return rothe_error_val, coeffs_new
        else:
            return rothe_error_val

    rothe_vg_jit = jax.jit(jax.value_and_grad(rothe_error, argnums=0), static_argnames=("SHH2", "return_cnew"))
    rothe_error_jit = jax.jit(rothe_error, static_argnames=("SHH2", "return_cnew"))
    return rothe_error_jit, rothe_vg_jit
