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
    regularization_lambda: float = 1e-12


# ---------------------------------------------------------------------------
# Tanh-bounded parameter transform
# ---------------------------------------------------------------------------


def compute_tanh_bounds(params_init_flat, params_shape, p=0.1, q=0.01, a_min=1e-4):
    """Compute element-wise (lb, ub) from the initial parameter vector.

    lb_j = init_j - p * |init_j| - q
    ub_j = init_j + p * |init_j| + q

    The lower bound for width parameters (index k=0 in the last axis of
    the (n, D, 4) layout) is clamped to *a_min* so that Gaussian widths
    stay strictly positive.
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
    return initial_err / 100  #


def maybe_apply_bb1_scaling(params_old, params_oldold, g_only_params, t, hess_inv_default):
    """
    Optionally apply BB1 (Barzilai-Borwein type 1) scaling to the
    initial inverse Hessian in parameter space.

    If the curvature condition is satisfied by the previous two steps,
    return BB1 * I. Otherwise, return the default.
    """
    if params_oldold is None:
        return hess_inv_default

    # Flatten p's
    params_old_flat = np.asarray(params_old).ravel()
    params_oldold_flat = np.asarray(params_oldold).ravel()

    # Step in parameter space.
    s_vec = params_old_flat - params_oldold_flat

    # Gradients at both previous solved steps in p-space (independent of start guess)
    grad_old_params_flat = np.asarray(g_only_params(params_old_flat))
    grad_oldold_params_flat = np.asarray(g_only_params(params_oldold_flat))

    # Gradient differences in parameter space.
    y_vec = grad_old_params_flat - grad_oldold_params_flat

    denom = float(s_vec @ y_vec)
    num = float(s_vec @ s_vec)
    BB_def = 1e8
    if not np.isfinite(denom) or not np.isfinite(num):
        BB1 = BB_def
    if denom <= 0.0:
        BB1 = BB_def
    else:
        denom_safe = denom
        BB1 = num / denom_safe
    if not np.isfinite(BB1) or BB1 <= 0.0:
        BB1 = BB_def

    BB1_max = BB_def
    if BB1 > BB1_max:
        print(f"{t} BB1={BB1:.3e} exceeds cap {BB1_max:.3e}, clamping.")
        BB1 = BB1_max

    dim = hess_inv_default.shape[0]
    return BB1 * np.eye(dim)
