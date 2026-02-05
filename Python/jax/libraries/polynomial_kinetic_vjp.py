"""
Custom VJP implementation for polynomial kinetic cross terms.

This provides analytical gradients for the polynomial kinetic cross terms,
avoiding autodiff through the moment calculation and index gathering.

NOTE on complex gradients in JAX:
JAX uses the convention where for f(z) with complex z:
- For product f = a * b: grad_a = g * b, grad_b = g * a  (NO conjugates)
- For einsum "ijm,m->ij": grad_poly[m] = sum_ij g[i,j] * term[i,j,m]  (NO conjugates)
- For hermitian h = VT + conj(VT).T: g_VT = g + conj(g).T

This is different from the Wirtinger derivative convention but matches JAX.
"""

from functools import partial

import jax
import jax.numpy as jnp

dtype_real = jnp.float64
dtype_complex = jnp.complex128


@partial(jax.custom_vjp, nondiff_argnums=(0, 1, 2, 3))
def compute_polynomial_kinetic_cross(
    dm,  # (D, M) int32 - nondiff
    idx,  # (D, M) int32 - nondiff
    idx1,  # (D, M) int32 - nondiff
    idx2,  # (D, M) int32 - nondiff
    moments,  # (n, n, D, K+1) complex
    tilde_k,  # (n, n, D) complex
    tilde_y,  # (n, n, D) complex
    beta_j,  # (n, n, D) complex
    S,  # (n, n) complex
    poly_vals,  # (M,) complex
):
    """
    Compute polynomial kinetic cross term VT + VT†.

    This is the forward computation with analytical gradient support.
    """
    n = moments.shape[0]
    D = moments.shape[2]
    M = poly_vals.shape[0]

    # Broadcast indices to (n, n, D, M)
    idx_full = jnp.broadcast_to(idx, (n, n, D, M))
    idx1_full = jnp.broadcast_to(idx1, (n, n, D, M))
    idx2_full = jnp.broadcast_to(idx2, (n, n, D, M))

    # Gather moments
    Mn = jnp.take_along_axis(moments, idx_full, axis=-1)  # (n,n,D,M)
    Mn1 = jnp.take_along_axis(moments, idx1_full, axis=-1)
    Mn2 = jnp.take_along_axis(moments, idx2_full, axis=-1)

    # Zero invalid entries
    dm_bc = dm[None, None, :, :]  # (1, 1, D, M)
    Mn1 = jnp.where(dm_bc >= 1, Mn1, 0.0)
    Mn2 = jnp.where(dm_bc >= 2, Mn2, 0.0)

    # Safe division
    safe_Mn = jnp.where(jnp.abs(Mn) > 0, Mn, 1.0)
    inv_Mn = jnp.where(jnp.abs(Mn) > 0, 1.0 / safe_Mn, 0.0)
    r1 = Mn1 * inv_Mn  # (n,n,D,M)
    r2 = Mn2 * inv_Mn  # (n,n,D,M)

    nvec = dm.astype(dtype_real)[None, None, :, :]  # (1,1,D,M)

    k = tilde_k[..., None]  # (n,n,D,1)
    y = tilde_y[..., None]  # (n,n,D,1)
    bj = beta_j[..., None]  # (n,n,D,1)

    base = -(2 * (k**2) * (y**2) - k)  # (n,n,D,1)
    vt_coeff = (
        base - 2 * nvec * k * y * bj * r1 - 0.5 * nvec * (nvec - 1) * (bj**2) * r2
    )  # (n,n,D,M)

    Mprod = jnp.prod(Mn, axis=2)  # (n,n,M)
    sum_coeff = jnp.sum(vt_coeff, axis=2)  # (n,n,M)

    term = S[..., None] * (Mprod * sum_coeff)  # (n,n,M)
    VT = jnp.einsum("ijm,m->ij", term, poly_vals)  # (n,n)

    return VT + VT.conj().T


def _polynomial_kinetic_fwd(
    dm, idx, idx1, idx2, moments, tilde_k, tilde_y, beta_j, S, poly_vals
):
    """Forward pass: compute result and save residuals for backward."""
    n = moments.shape[0]
    D = moments.shape[2]
    M = poly_vals.shape[0]

    # Broadcast indices
    idx_full = jnp.broadcast_to(idx, (n, n, D, M))
    idx1_full = jnp.broadcast_to(idx1, (n, n, D, M))
    idx2_full = jnp.broadcast_to(idx2, (n, n, D, M))

    # Gather moments
    Mn = jnp.take_along_axis(moments, idx_full, axis=-1)
    Mn1 = jnp.take_along_axis(moments, idx1_full, axis=-1)
    Mn2 = jnp.take_along_axis(moments, idx2_full, axis=-1)

    dm_bc = dm[None, None, :, :]
    Mn1 = jnp.where(dm_bc >= 1, Mn1, 0.0)
    Mn2 = jnp.where(dm_bc >= 2, Mn2, 0.0)

    safe_Mn = jnp.where(jnp.abs(Mn) > 0, Mn, 1.0)
    inv_Mn = jnp.where(jnp.abs(Mn) > 0, 1.0 / safe_Mn, 0.0)
    r1 = Mn1 * inv_Mn
    r2 = Mn2 * inv_Mn

    nvec = dm.astype(dtype_real)[None, None, :, :]

    k = tilde_k[..., None]
    y = tilde_y[..., None]
    bj = beta_j[..., None]

    base = -(2 * (k**2) * (y**2) - k)
    vt_coeff = base - 2 * nvec * k * y * bj * r1 - 0.5 * nvec * (nvec - 1) * (bj**2) * r2

    Mprod = jnp.prod(Mn, axis=2)
    sum_coeff = jnp.sum(vt_coeff, axis=2)

    term = S[..., None] * (Mprod * sum_coeff)
    VT = jnp.einsum("ijm,m->ij", term, poly_vals)
    result = VT + VT.conj().T

    # Save residuals for backward pass
    residuals = (
        moments,
        tilde_k,
        tilde_y,
        beta_j,
        S,
        poly_vals,
        Mn,
        Mn1,
        Mn2,
        r1,
        r2,
        inv_Mn,
        Mprod,
        sum_coeff,
        term,
        VT,
    )

    return result, residuals


def _polynomial_kinetic_bwd(dm, idx, idx1, idx2, residuals, g):
    """
    Backward pass: compute VJP analytically.

    JAX complex VJP rules (verified empirically):
    - For product a*b: grad_a = g * b, grad_b = g * a  (NO conjugates)
    - For einsum: grad = einsum with same structure (NO conjugates on data)
    - For conj(x): grad = conj(g)
    - For x.T: grad = g.T

    g: cotangent w.r.t. output (n, n) complex

    Returns: tuple of cotangents for (moments, tilde_k, tilde_y, beta_j, S, poly_vals)
    """
    (
        moments,
        tilde_k,
        tilde_y,
        beta_j,
        S,
        poly_vals,
        Mn,
        Mn1,
        Mn2,
        r1,
        r2,
        inv_Mn,
        Mprod,
        sum_coeff,
        term,
        VT,
    ) = residuals

    n = tilde_k.shape[0]
    D = tilde_k.shape[2]
    M = poly_vals.shape[0]

    nvec = dm.astype(dtype_real)[None, None, :, :]  # (1,1,D,M)
    dm_bc = dm[None, None, :, :]

    k = tilde_k[..., None]  # (n,n,D,1)
    y = tilde_y[..., None]
    bj = beta_j[..., None]

    # Result = VT + conj(VT).T
    # For h(VT) = VT + conj(VT).T:
    #   ∂h/∂VT[i,j] contributes to Result[i,j] with weight 1
    #   ∂h/∂VT*[i,j] = conj of VT at (j,i) contributes to Result[j,i]
    # VJP: g_VT[i,j] = g[i,j] (from VT term) + conj(g[j,i]) (from conj(VT.T) term)
    g_VT = g + jnp.conj(g.T)  # (n, n)

    # VT = einsum("ijm,m->ij", term, poly_vals)
    # JAX VJP (empirically verified):
    # grad_poly_vals = einsum("ij,ijm->m", g_VT, term)  (NO conjugates)
    # grad_term = g_VT[..., None] * poly_vals (for each m)? No...
    # Actually: grad_term[i,j,m] = g_VT[i,j] * poly_vals[m] (NO conjugates)
    grad_poly_vals = jnp.einsum("ij,ijm->m", g_VT, term)
    g_term = g_VT[..., None] * poly_vals[None, None, :]  # (n,n,M)

    # term = S[..., None] * (Mprod * sum_coeff)
    # For product a*b: grad_a = g*b, grad_b = g*a
    # Here: term = S * (Mprod * sum_coeff)
    # grad_S[i,j] = sum_m g_term[i,j,m] * (Mprod * sum_coeff)[i,j,m]
    grad_S = jnp.sum(g_term * (Mprod * sum_coeff), axis=-1)  # (n,n)

    # g_(Mprod * sum_coeff) = g_term * S[..., None]
    g_MS = g_term * S[..., None]  # (n,n,M)

    # Mprod * sum_coeff -> product rule (no conjugates)
    # grad_Mprod = g_MS * sum_coeff
    # grad_sum_coeff = g_MS * Mprod
    g_Mprod = g_MS * sum_coeff  # (n,n,M)
    g_sum_coeff = g_MS * Mprod  # (n,n,M)

    # sum_coeff = sum_d vt_coeff over axis=2
    # grad_vt_coeff[i,j,d,m] = g_sum_coeff[i,j,m] (broadcast over d)
    g_vt_coeff = g_sum_coeff[:, :, None, :]  # (n,n,1,M) -> broadcasts to (n,n,D,M)

    # vt_coeff = base - 2*nvec*k*y*bj*r1 - 0.5*nvec*(nvec-1)*bj^2*r2
    # where base = -(2*k^2*y^2 - k) = k - 2*k^2*y^2

    # grad_base (summed over M dimension since base is (n,n,D,1))
    g_base_summed = jnp.sum(g_vt_coeff, axis=-1, keepdims=True)  # (n,n,D,1)

    # base = k - 2*k^2*y^2
    # ∂base/∂k = 1 - 4*k*y^2
    # ∂base/∂y = -4*k^2*y
    g_k_from_base = g_base_summed * (1 - 4 * k * (y**2))  # (n,n,D,1)
    g_y_from_base = g_base_summed * (-4 * (k**2) * y)  # (n,n,D,1)

    # Term: -2*nvec*k*y*bj*r1
    # ∂/∂k = -2*nvec*y*bj*r1
    # ∂/∂y = -2*nvec*k*bj*r1
    # ∂/∂bj = -2*nvec*k*y*r1
    # ∂/∂r1 = -2*nvec*k*y*bj
    g_k_from_r1term = g_vt_coeff * (-2 * nvec * y * bj * r1)  # (n,n,D,M)
    g_y_from_r1term = g_vt_coeff * (-2 * nvec * k * bj * r1)
    g_bj_from_r1term = g_vt_coeff * (-2 * nvec * k * y * r1)
    g_r1 = g_vt_coeff * (-2 * nvec * k * y * bj)  # (n,n,D,M)

    # Term: -0.5*nvec*(nvec-1)*bj^2*r2
    # ∂/∂bj = -nvec*(nvec-1)*bj*r2
    # ∂/∂r2 = -0.5*nvec*(nvec-1)*bj^2
    g_bj_from_r2term = g_vt_coeff * (-nvec * (nvec - 1) * bj * r2)  # (n,n,D,M)
    g_r2 = g_vt_coeff * (-0.5 * nvec * (nvec - 1) * (bj**2))  # (n,n,D,M)

    # Sum up k, y, bj gradients
    grad_tilde_k = jnp.sum(g_k_from_base, axis=-1) + jnp.sum(
        g_k_from_r1term, axis=-1
    )  # (n,n,D)
    grad_tilde_y = jnp.sum(g_y_from_base, axis=-1) + jnp.sum(
        g_y_from_r1term, axis=-1
    )  # (n,n,D)
    grad_beta_j = jnp.sum(g_bj_from_r1term, axis=-1) + jnp.sum(
        g_bj_from_r2term, axis=-1
    )  # (n,n,D)

    # Now backprop through r1 = Mn1 * inv_Mn, r2 = Mn2 * inv_Mn
    # r1 = Mn1 / Mn (where Mn != 0)
    # For quotient a/b: grad_a = g/b, grad_b = -g*a/b^2 = -g*r/b
    # r1 = Mn1 * inv_Mn
    # grad_Mn1 = g_r1 * inv_Mn  (where dm >= 1)
    # grad_Mn from r1: g_r1 * Mn1 * (-inv_Mn^2) = -g_r1 * r1 * inv_Mn
    g_Mn1 = jnp.where(dm_bc >= 1, g_r1 * inv_Mn, 0.0)
    g_Mn_from_r1 = jnp.where(jnp.abs(Mn) > 0, -g_r1 * r1 * inv_Mn, 0.0)

    g_Mn2 = jnp.where(dm_bc >= 2, g_r2 * inv_Mn, 0.0)
    g_Mn_from_r2 = jnp.where(jnp.abs(Mn) > 0, -g_r2 * r2 * inv_Mn, 0.0)

    # Mprod = prod_d Mn[..., d, :]
    # ∂Mprod/∂Mn[d] = Mprod / Mn[d] (where Mn[d] != 0)
    safe_Mn_prod = jnp.where(jnp.abs(Mn) > 0, Mn, 1.0)
    inv_Mn_prod = jnp.where(jnp.abs(Mn) > 0, 1.0 / safe_Mn_prod, 0.0)
    # g_Mprod is (n,n,M), need to broadcast to (n,n,D,M)
    # grad_Mn[i,j,d,m] = g_Mprod[i,j,m] * Mprod[i,j,m] / Mn[i,j,d,m]
    g_Mn_from_Mprod = g_Mprod[:, :, None, :] * Mprod[:, :, None, :] * inv_Mn_prod  # (n,n,D,M)

    # Total gradient for Mn
    g_Mn = g_Mn_from_r1 + g_Mn_from_r2 + g_Mn_from_Mprod  # (n,n,D,M)

    # Now scatter g_Mn, g_Mn1, g_Mn2 back to moments
    # moments is (n, n, D, K+1), indices are (D, M)
    K_plus_1 = moments.shape[-1]

    # Create one-hot encoding: (D, M, K+1) -> which moment each (d,m) maps to
    one_hot_idx = jax.nn.one_hot(idx, K_plus_1, dtype=dtype_complex)  # (D, M, K+1)
    one_hot_idx1 = jax.nn.one_hot(idx1, K_plus_1, dtype=dtype_complex)
    one_hot_idx2 = jax.nn.one_hot(idx2, K_plus_1, dtype=dtype_complex)

    # g_Mn is (n, n, D, M), one_hot is (D, M, K+1)
    # We want: grad_moments[i, j, d, k] = sum_m g_Mn[i, j, d, m] * one_hot[d, m, k]
    grad_moments_from_Mn = jnp.einsum("ijdm,dmk->ijdk", g_Mn, one_hot_idx)
    grad_moments_from_Mn1 = jnp.einsum("ijdm,dmk->ijdk", g_Mn1, one_hot_idx1)
    grad_moments_from_Mn2 = jnp.einsum("ijdm,dmk->ijdk", g_Mn2, one_hot_idx2)

    grad_moments = grad_moments_from_Mn + grad_moments_from_Mn1 + grad_moments_from_Mn2

    return (grad_moments, grad_tilde_k, grad_tilde_y, grad_beta_j, grad_S, grad_poly_vals)


# ============================================================================
# Custom VJP Registration (currently disabled - using autodiff instead)
# ============================================================================
# The custom VJP is correct but currently slower than JAX autodiff for small
# problem sizes. To enable, uncomment the line below:
#
# compute_polynomial_kinetic_cross.defvjp(_polynomial_kinetic_fwd, _polynomial_kinetic_bwd)
#
# Note: For real speedups, a full-chain VJP from params -> output is needed,
# similar to calculate_Gaussian_expectation_values in gaussian_potential_helpers.py
# ============================================================================
