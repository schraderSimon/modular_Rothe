"""
Optimization helper functions for the Rothe solver.

Contains utility functions for parameter transformation, objective construction,
Hessian approximation, and BFGS initialization strategies.
"""

import numpy as np

import jax
import jax.numpy as jnp

from .calculate_Hessian_coefficients import calculate_Hessian_matrix
from .general_wf import ND_potentials


@jax.jit
def evaluate_function_1000_times(func, x):
    import time

    start = time.time()
    for _ in range(1000):
        func(x)
    end = time.time()
    time_taken = end - start
    print(f"Time taken for 1000 evaluations: {time_taken:.2f} seconds")
    return time_taken


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
    scale_arr = np.einsum(
        "...k,k->...k", np.ones_like(grad_p), scale_types
    )  # last axis broadcasts
    S_flat = scale_arr.ravel()
    S_flat_jax = jnp.asarray(S_flat)
    return S_flat, S_flat_jax


def params_flat_to_y_flat(params_flat, S_flat):
    """
    Change of basis params = S y  ⇒  y = S^{-1} params.

    params_flat: np.ndarray, flattened params
    S_flat: np.ndarray, same shape, diagonal scaling
    """
    return params_flat / S_flat


def y_flat_to_params_flat(y_flat, S_flat):
    """
    Inverse change of basis y → params = S y.

    y_flat: np.ndarray, flattened y
    S_flat: np.ndarray, same shape, diagonal scaling
    """
    return y_flat * S_flat


def make_objective_params(
    rothe_error_and_gradient, params_shape, params_old, coeffs_old, SHH2, t, dt
):
    """
    Wrap rothe_error_and_gradient into a flat params-space objective.

    Returns:
        f_and_g_params(theta_flat) -> (float value, np.ndarray grad_flat)
        g_only_params(theta_flat)  -> grad_flat
    """

    def f_and_g_params(theta_flat):
        params_new = jnp.asarray(theta_flat).reshape(params_shape)
        val, g = rothe_error_and_gradient(params_new, params_old, coeffs_old, SHH2, t, dt)
        return float(val), np.asarray(jnp.ravel(g))

    def g_only_params(theta_flat):
        _, g = f_and_g_params(theta_flat)
        return g

    return f_and_g_params, g_only_params


def make_objective_y(
    rothe_error_and_gradient,
    params_shape,
    params_old,
    coeffs_old,
    SHH2,
    t,
    dt,
    S_flat_jax,
):
    """
    Wrap rothe_error_and_gradient into y-space objective, with params = S y.

    Returns:
        f_and_g_y(theta_y_flat) -> (float value, np.ndarray grad_y_flat)
    """

    def f_and_g_y(theta_y_flat):
        theta_y = jnp.asarray(theta_y_flat)
        params_flat = theta_y * S_flat_jax  # params = S * y
        params_new = params_flat.reshape(params_shape)
        val, grad_params = rothe_error_and_gradient(
            params_new, params_old, coeffs_old, SHH2, t, dt
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


def maybe_apply_bb1_scaling(
    params_old,
    params_oldold,
    T,
    T_inv,
    grad_theta0_flat,
    g_only_params,
    t,
    hess_inv_default,
):

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

    print(f"{t} fulfills curvature condition (theta-space): BB1={BB1}")
    dim = hess_inv_default.shape[0]
    return BB1 * np.eye(dim)


def get_approx_Hessian(params_init, coeffs):
    wf_obj = ND_potentials(params_init, 4)
    S = wf_obj.calculate_S()
    moments = wf_obj.moments
    hess_new = calculate_Hessian_matrix(params_init, S, moments, coeffs)
    n, D, _ = params_init.shape
    h_p_6d = jnp.transpose(hess_new, (0, 2, 4, 1, 3, 5))  # (n, D, 4, n, D, 4)
    i = np.arange(n)
    mask = i[:, None] == i[None, :]  # shape (n, n)

    # reshape to broadcast over (n, D, 4, n, D, 4)
    mask = mask[:, None, None, :, None, None]

    h_p_6d *= mask
    H_p = h_p_6d.reshape(n * D * 4, n * D * 4)
    # H_p = H_p + H_p.T  # Ensure Hermiticity

    return H_p


def hessian_cholesky(hessian):
    """Return a factor L such that H ≈ L^T L.

    Note:
        This currently clips eigenvalues to enforce positive definiteness.
        The function is conservative about numerical stability.

    Caveat:
        The dense-metric path in `RotheSolver.find_next_timestep_solution` is
        experimental; if you rely on this for preconditioning, confirm the
        scaling is actually enabled.
    """
    hessian_np = np.asarray(hessian, dtype=np.float64)
    w, Q = np.linalg.eigh(hessian_np)
    # Clip eigenvalues to avoid insane conditioning / negatives
    lam_med = np.median(w[w > 0])  # median of positive ones
    lam_min = lam_med * 1e-3
    w_clipped = np.clip(w, lam_min, None)
    # w_clipped = np.ones_like(w)  # TEMPORARY: use identity
    # Build L such that H_p ≈ L^T L
    sqrt_w = np.sqrt(w_clipped)
    L = sqrt_w[None, :] * Q.T
    return np.eye(L.shape[0])  # TEMPORARY: use identity


def make_objective_z(
    rothe_error_and_gradient,
    params_shape,
    params_old,
    coeffs_old,
    SHH2,
    t,
    dt,
    L_inv,
):
    """
    z-space objective: params = L^{-1} z, g_z = L^{-T} grad_params
    """
    L_inv_T = L_inv.T

    def f_and_g_z(z_flat):
        z = jnp.asarray(z_flat)
        # p = L^{-1} z
        params_flat = L_inv @ np.asarray(z)  # use np matmul, that's fine here
        params_new = params_flat.reshape(params_shape)
        val, grad_params = rothe_error_and_gradient(
            params_new, params_old, coeffs_old, SHH2, t, dt
        )
        grad_params_flat = np.asarray(jnp.ravel(grad_params))
        # g_z = L^{-T} g_p
        g_z_flat = L_inv_T @ grad_params_flat
        return float(val), g_z_flat

    return f_and_g_z
