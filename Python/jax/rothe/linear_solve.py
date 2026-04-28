"""Shared linear solve and scoring helpers for Rothe coefficient updates."""

import numpy as np

import jax.numpy as jnp
from jax.scipy.sparse.linalg import cg as jax_cg


def solve_and_score_jax(A_dagger_A, rho_mat, B_dagger_B, coeffs_old, lambda_):
    """Solve coefficients and compute Rothe error with JAX arrays."""
    rho_vec = jnp.conj(rho_mat).T @ coeffs_old
    A_reg = A_dagger_A + jnp.eye(A_dagger_A.shape[0], dtype=A_dagger_A.dtype) * lambda_
    # Use iterative CG to avoid dense cuSolver calls that can fail to
    # initialize on some GPU setups.
    diag = jnp.diag(A_reg)
    abs_diag = jnp.abs(diag)
    if jnp.issubdtype(A_reg.dtype, jnp.complexfloating):
        one = jnp.array(1.0 + 0.0j, dtype=A_reg.dtype)
    else:
        one = jnp.array(1.0, dtype=A_reg.dtype)
    safe_diag = jnp.where(abs_diag > 1e-14, diag, one)

    n = int(A_reg.shape[0])
    coeffs_new, _ = jax_cg(
        lambda x: A_reg @ x, rho_vec, M=lambda x: x / safe_diag, tol=1e-12, atol=1e-15, maxiter=max(50, 10 * n)
    )
    rothe_error_val = jnp.real(jnp.conj(coeffs_old) @ B_dagger_B @ coeffs_old - jnp.conj(rho_vec).T @ coeffs_new)
    return rothe_error_val, coeffs_new


def solve_and_score_numpy(A_dagger_A, rho_mat, B_dagger_B, coeffs_old, lambda_):
    """Solve coefficients and compute Rothe error with NumPy arrays."""
    A_np = np.asarray(A_dagger_A)
    rho_np = np.asarray(rho_mat)
    B_np = np.asarray(B_dagger_B)
    coeffs_old_np = np.asarray(coeffs_old)

    rho_vec = np.conj(rho_np).T @ coeffs_old_np
    A_reg = A_np + np.eye(A_np.shape[0], dtype=A_np.dtype) * lambda_
    coeffs_new = np.linalg.solve(A_reg, rho_vec)
    rothe_error_val = np.real(np.conj(coeffs_old_np) @ B_np @ coeffs_old_np - np.conj(rho_vec).T @ coeffs_new)
    return float(rothe_error_val), coeffs_new
