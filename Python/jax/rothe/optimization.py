"""Optimization helper functions for the Rothe solver."""

from dataclasses import dataclass
from typing import Any

import numpy as np

import jax.numpy as jnp


@dataclass(frozen=True)
class RotheObjectiveContext:
    """Shared inputs passed together across Rothe objective evaluations."""

    params_old: Any
    coeffs_old: Any
    SHH2: Any
    t: Any
    dt: Any
    params_frozen: Any = None
    SHH2_oldold: Any = None
    regularization_lambda: float | None = None


# ---------------------------------------------------------------------------
# Tanh-bounded parameter transform
# ---------------------------------------------------------------------------


def compute_tanh_bounds(params_init_flat, params_shape, p=0.1, q=0.01, a_min=1e-2, a_max=1e2):
    """Compute element-wise (lb, ub) from the initial parameter vector.

    lb_j = init_j - p * |init_j| - q
    ub_j = init_j + p * |init_j| + q

    Bounds are symmetric around the initial point for all coordinates.
    Additionally, the real width part ``a`` (k=0 in ``(n, D, 4)``) is clipped
    to ``[a_min, a_max]`` so Gaussian widths stay in a stable positive range.
    """
    params_init_flat = np.asarray(params_init_flat, dtype=float)
    half_width = p * np.abs(params_init_flat) + q

    lb = params_init_flat - half_width
    ub = params_init_flat + half_width

    # Clamp lower bound of width parameter a (k=0) to a_min.
    a_mask = np.zeros(params_shape, dtype=bool)
    a_mask[:, :, 0] = True
    a_mask_flat = a_mask.ravel()
    lb[a_mask_flat] = np.maximum(lb[a_mask_flat], a_min)
    ub[a_mask_flat] = np.minimum(ub[a_mask_flat], a_max)
    return lb, ub


def tanh_forward(alpha_prime, lb, ub):
    """Map unconstrained α' → bounded α ∈ (lb, ub).

    α = lb + (tanh(α') + 1) * (ub - lb) / 2
    """
    return lb + (np.tanh(alpha_prime) + 1.0) * (ub - lb) / 2.0


def tanh_inverse(alpha, lb, ub):
    """Map bounded α → unconstrained α'.

    α' = atanh(2*(α - lb) / (ub - lb) - 1)

    At the initial point (midpoint of [lb, ub]) this returns 0.
    """
    z = 2.0 * (alpha - lb) / (ub - lb) - 1.0
    z = np.clip(z, -1 + 1e-12, 1 - 1e-12)
    return np.arctanh(z)


def tanh_forward_jacobian(alpha_prime, lb, ub):
    """Element-wise dα/dα' = (ub - lb)/2 · sech²(α')."""
    return (ub - lb) / 2.0 * (1.0 - np.tanh(alpha_prime) ** 2)


def make_objective_params(rothe_error_and_gradient, params_shape, context):
    """
    Wrap rothe_error_and_gradient into a flat params-space objective.

    The required ``context`` carries the shared objective inputs
    (params_old, coeffs_old, SHH2, t, dt, params_frozen, SHH2_oldold).

    Returns:
        f_and_g_params(theta_flat) -> (float value, np.ndarray grad_flat)
        g_only_params(theta_flat)  -> grad_flat
    """
    params_old = context.params_old
    coeffs_old = context.coeffs_old
    SHH2 = context.SHH2
    t = context.t
    dt = context.dt
    params_frozen = context.params_frozen
    SHH2_oldold = context.SHH2_oldold
    regularization_lambda = context.regularization_lambda
    if regularization_lambda is None:
        raise ValueError("regularization_lambda must be provided in RotheObjectiveContext.")

    def f_and_g_params(theta_flat):
        params_new = jnp.asarray(theta_flat).reshape(params_shape)
        val, g = rothe_error_and_gradient(
            params_new,
            params_old,
            coeffs_old,
            SHH2,
            t,
            dt,
            params_frozen,
            lambda_=regularization_lambda,
            SHH2_oldold=SHH2_oldold,
        )
        val_float = float(val)
        grad_flat = np.asarray(jnp.ravel(g))
        if not np.isfinite(val_float) or not np.all(np.isfinite(grad_flat)):
            return 1e30, np.zeros_like(np.asarray(theta_flat), dtype=float)
        return val_float, grad_flat

    def g_only_params(theta_flat):
        _, g = f_and_g_params(theta_flat)
        return g

    return f_and_g_params, g_only_params


def make_objective_tanh(f_and_g_params, lb, ub):
    """Wrap a flat params-space objective into the tanh-bounded space.

    The optimizer works with unconstrained α'; internally this applies
    α = tanh_forward(α', lb, ub) and chains the gradient via the Jacobian.

    Returns:
        f_and_g_tanh(alpha_prime_flat) -> (float value, np.ndarray grad_flat)
    """
    lb = np.asarray(lb, dtype=float)
    ub = np.asarray(ub, dtype=float)

    def f_and_g_tanh(alpha_prime_flat):
        alpha_prime = np.asarray(alpha_prime_flat, dtype=float)
        alpha = tanh_forward(alpha_prime, lb, ub)
        val, grad_alpha = f_and_g_params(alpha)
        jac = tanh_forward_jacobian(alpha_prime, lb, ub)
        grad_prime = grad_alpha * jac
        return val, grad_prime

    return f_and_g_tanh


def compute_gtol(initial_err):
    """
    Stopping criterion for BFGS.
    """
    return initial_err / 10  #


def make_scaled_identity(scale, dim):
    """Return ``scale * I`` with shape ``(dim, dim)``."""
    return float(scale) * np.eye(int(dim), dtype=float)


def compute_bb1_inverse_hessian_scale(
    params_old, params_oldold, g_only_params, *, default_scale=1.0, max_scale=1e8, t=None
):
    """Compute BB1 scalar for inverse-Hessian initialization.

    BB1 uses
        s = x_k - x_{k-1}
        y = g_k - g_{k-1}
        H0 = (s^T s / s^T y) I

    If BB1 is not usable (missing previous step, bad curvature, non-finite
    values), ``default_scale`` is returned.
    """
    default_scale = float(default_scale)
    if not np.isfinite(default_scale) or default_scale <= 0.0:
        default_scale = 1.0

    if params_oldold is None:
        return default_scale

    params_old_flat = np.asarray(params_old).ravel()
    params_oldold_flat = np.asarray(params_oldold).ravel()
    s_vec = params_old_flat - params_oldold_flat

    grad_old_params_flat = np.asarray(g_only_params(params_old_flat))
    grad_oldold_params_flat = np.asarray(g_only_params(params_oldold_flat))
    y_vec = grad_old_params_flat - grad_oldold_params_flat

    denom = float(s_vec @ y_vec)
    num = float(s_vec @ s_vec)
    if not np.isfinite(denom) or not np.isfinite(num) or num <= 0.0 or denom <= 0.0:
        return default_scale

    bb1_scale = num / denom
    if not np.isfinite(bb1_scale) or bb1_scale <= 0.0:
        return default_scale

    max_scale = float(max_scale)
    if np.isfinite(max_scale) and max_scale > 0.0 and bb1_scale > max_scale:
        if t is not None:
            print(f"{t} BB1={bb1_scale:.3e} exceeds cap {max_scale:.3e}, clamping.")
        bb1_scale = max_scale

    return bb1_scale
