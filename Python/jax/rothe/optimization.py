"""
Optimization helper functions for the Rothe solver.

Contains utility functions for parameter transformation, objective construction,
and BFGS initialization strategies (including BB1 scaling).
"""

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


def build_diagonal_scaling_from_grad(grad_p, epsilon=1e-12):
    """
    Build diagonal scaling S from p-space gradient.

    grad_p: np.ndarray, shape (n, D, 4)
        Gradient in p-space reshaped like p.
    Returns:
        S_flat:     np.ndarray, shape (n*D*4,)   (diagonal of S)
        S_flat_jax: jnp.ndarray, same data, for JAX ops
    """
    abs_grad = np.abs(grad_p)
    # RMS per parameter type (last axis: a, b, mu, p)
    group_rms = np.sqrt(np.mean(abs_grad**2, axis=(0, 1)))  # shape (4,)

    # Target RMS: average of the four
    desired_RMS = np.mean(group_rms)

    # Scale factors: s_k = avg / RMS_k, avoid division by zero
    scale_types = np.where(group_rms > epsilon, desired_RMS / group_rms, 1.0)  # shape (4,)

    # Broadcast to full shape (n, D, 4)
    scale_arr = np.einsum("...k,k->...k", np.ones_like(grad_p), scale_types)  # last axis broadcasts
    S_flat = scale_arr.ravel()
    S_flat_jax = jnp.asarray(S_flat)
    return S_flat, S_flat_jax


def params_flat_to_y_flat(params_flat, S_flat):
    """
    Change of basis params = S y  =>  y = S^{-1} params.
    """
    return params_flat / S_flat


def y_flat_to_params_flat(y_flat, S_flat):
    """
    Inverse change of basis y -> params = S y.
    """
    return y_flat * S_flat


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

    def f_and_g_params(theta_flat):
        params_new = jnp.asarray(theta_flat).reshape(params_shape)
        val, g = rothe_error_and_gradient(
            params_new, params_old, coeffs_old, SHH2, t, dt, params_frozen, SHH2_oldold=SHH2_oldold
        )
        return float(val), np.asarray(jnp.ravel(g))

    def g_only_params(theta_flat):
        _, g = f_and_g_params(theta_flat)
        return g

    return f_and_g_params, g_only_params


def make_objective_y(rothe_error_and_gradient, params_shape, S_flat_jax, context):
    """
    Wrap rothe_error_and_gradient into y-space objective, with params = S y.

    The required ``context`` carries the shared objective inputs
    (params_old, coeffs_old, SHH2, t, dt, params_frozen, SHH2_oldold).

    Returns:
        f_and_g_y(theta_y_flat) -> (float value, np.ndarray grad_y_flat)
    """
    params_old = context.params_old
    coeffs_old = context.coeffs_old
    SHH2 = context.SHH2
    t = context.t
    dt = context.dt
    params_frozen = context.params_frozen
    SHH2_oldold = context.SHH2_oldold

    def f_and_g_y(theta_y_flat):
        theta_y = jnp.asarray(theta_y_flat)
        params_flat = theta_y * S_flat_jax  # params = S * y
        params_new = params_flat.reshape(params_shape)
        val, grad_params = rothe_error_and_gradient(
            params_new, params_old, coeffs_old, SHH2, t, dt, params_frozen, SHH2_oldold=SHH2_oldold
        )
        grad_params_flat = jnp.ravel(grad_params)
        g_y_flat = grad_params_flat * S_flat_jax  # g_y = S * grad_params
        return float(val), np.asarray(g_y_flat)

    return f_and_g_y


def compute_gtol(initial_err):
    """
    Stopping criterion for BFGS.
    """
    return max(initial_err / 10.0, 1e-10)


def make_initial_hessian_y(dim_y, grad_y0_norm):
    """
    Simple scaled identity initial inverse Hessian in y-space.
    """
    return np.eye(dim_y) / grad_y0_norm


def maybe_apply_bb1_scaling(params_old, params_oldold, T, T_inv, grad_theta0_flat, g_only_params, t, hess_inv_default):
    """
    Optionally apply BB1 (Barzilai-Borwein type 1) scaling to the
    initial inverse Hessian in theta-space.

    If the curvature condition is satisfied by the previous two steps,
    return BB1 * I. Otherwise, return the default.
    """
    if params_oldold is None:
        return hess_inv_default

    # Flatten p's
    params_old_flat = np.asarray(params_old).ravel()
    params_oldold_flat = np.asarray(params_oldold).ravel()

    # Step in p-space
    delta_params = params_old_flat - params_oldold_flat  # s_params

    # Step in theta-space: s_theta = T^{-1} delta_p
    s_theta = T_inv @ delta_params

    # Gradient at old-old point in p-space
    grad_oldold_params_flat = np.asarray(g_only_params(params_oldold_flat))

    # Gradient in theta at old-old: g_theta_oldold = T^T g_p_oldold
    g_theta_oldold = T.T @ grad_oldold_params_flat

    # Difference in theta-gradients
    grad_theta0_flat = np.asarray(grad_theta0_flat)
    y_vec = grad_theta0_flat - g_theta_oldold

    denom = float(s_theta @ y_vec)
    num = float(s_theta @ s_theta)

    if not np.isfinite(denom) or denom <= 0.0:
        print(f"{t} does not fulfill curvature condition: denom={denom}, using default.")
        return hess_inv_default

    BB1 = num / denom
    if not np.isfinite(BB1) or BB1 <= 0.0:
        print(f"{t} invalid BB1={BB1}, using default.")
        return hess_inv_default

    dim = hess_inv_default.shape[0]
    return BB1 * np.eye(dim)
