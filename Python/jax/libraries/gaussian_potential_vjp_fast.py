"""
Fast vectorized VJP for Gaussian potential expectation values.

This module provides two implementations:

1. calculate_Gaussian_expectation_values_fast:
   - Fully vectorized over all m potential terms
   - Best performance for moderate problem sizes
   - Memory: O(m * n² * D * 16 bytes) per tensor, ~15 tensors in backward pass

2. calculate_Gaussian_expectation_values_batched:
   - Batched version that processes m terms in chunks
   - Configurable memory limit via memory_limit_mb parameter
   - Falls back to fully vectorized when batch_size >= m
   - Speed vs memory trade-off: fewer batches = faster but more memory

Memory estimate:
    mem_gb ≈ 15 * m * n² * D * 16 / 1e9

Performance notes:
- For small/moderate m: Use fast version (single batch)
- For large m where memory is a concern: Use batched with appropriate limit
"""

from functools import partial

import jax
import jax.numpy as jnp

from .helpers import get_individual_params

dtype_complex = jnp.complex128
dtype_real = jnp.float64


# =============================================================================
# CORE HELPER FUNCTIONS
# =============================================================================


def get_gaussian_pair_terms(params):
    """
    Compute pairwise Gaussian terms for basis overlap.

    Args:
        params: (n, D, 4) basis parameters [a, b, mu, p]

    Returns:
        gauss_A, gauss_B, gauss_C: (n, n, D) tensors defining the pairwise Gaussian
    """
    a, b, mu, p = get_individual_params(params)
    alpha = a**2 + 1j * b  # (n, D)
    ai, aj = alpha.conj()[:, None, :], alpha[None, :, :]

    gauss_A = ai + aj  # (n, n, D)
    gauss_B = 2 * (ai * mu[:, None, :] + aj * mu[None, :, :]) + 1j * (
        p[None, :, :] - p[:, None, :]
    )
    gauss_C = -(ai * mu[:, None, :] ** 2 + aj * mu[None, :, :] ** 2) + 1j * (
        mu[:, None, :] * p[:, None, :] - mu[None, :, :] * p[None, :, :]
    )
    return gauss_A, gauss_B, gauss_C


def get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C):
    """Compute normalization terms from basis overlap."""
    A_diag = gauss_A.diagonal(axis1=0, axis2=1)  # (D, n)
    B_diag = gauss_B.diagonal(axis1=0, axis2=1)
    C_diag = gauss_C.diagonal(axis1=0, axis2=1)

    center_diag = B_diag / (2 * A_diag)
    constant_diag = C_diag + center_diag**2 * A_diag
    exponent_subtraction = (constant_diag - 0.5 * jnp.log(A_diag)).T  # (n, D)

    return exponent_subtraction[:, None, :] + exponent_subtraction[None, :, :]  # (n, n, D)


def combine_with_potential(gauss_A, gauss_B, gauss_C, alphaK, muK, pK):
    """
    Combine basis pair terms with potential Gaussian parameters.

    Adds potential Gaussian (alphaK, muK, pK) to basis pair (gauss_A, gauss_B, gauss_C).
    Works for any leading batch dimensions on alphaK, muK, pK.

    Args:
        gauss_A, gauss_B, gauss_C: (n, n, D) basis pair terms
        alphaK, muK, pK: (..., D) potential parameters (can have batch dims)

    Returns:
        allA, allB, allC: (..., n, n, D) combined terms
    """
    # Expand potential params for broadcasting with (n, n, D)
    expand = lambda x: x[..., None, None, :]  # (..., D) -> (..., 1, 1, D)
    aK, mK, pK = expand(alphaK), expand(muK), expand(pK)

    allA = gauss_A + aK
    allB = gauss_B + 2 * aK * mK + 1j * pK
    allC = gauss_C - aK * mK**2 - 1j * mK * pK
    return allA, allB, allC


def compute_contributions(allA, allB, allC, exp_substr_terms, linear_params):
    """
    Compute weighted Gaussian contributions.

    Args:
        allA, allB, allC: (m, n, n, D) combined Gaussian terms
        exp_substr_terms: (n, n, D) normalization terms
        linear_params: (m,) linear coefficients

    Returns:
        contribs: (m, n, n) weighted contributions per potential term
        centers: (m, n, n, D) Gaussian centers (needed for gradients)
    """
    centers = allB / (2 * allA)
    constants = allC + centers**2 * allA
    sum_exponents = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms[None, :, :, :]
    sum_exponents_ij = jnp.sum(sum_exponents, axis=-1)  # (m, n, n)

    contribs = jnp.exp(sum_exponents_ij) * linear_params[:, None, None]
    return contribs, centers


def compute_param_gradients(weighted, centers, allA, alpha, mu, p):
    """
    Compute gradients w.r.t. basis parameters from weighted contributions.

    Optimization: Contract over m first to get (n,n,D), then sum rows/cols.
    This creates smaller intermediates than direct (m,n,n,D) -> (n,D) contraction.

    Args:
        weighted: (m, n, n) weighted output gradients
        centers: (m, n, n, D) Gaussian centers
        allA: (m, n, n, D) combined A terms
        alpha, mu, p: (n, D) basis parameters

    Returns:
        Tuple: (grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_sum)
    """
    # Core quantities
    inv_A = 1.0 / allA
    centers_sq = centers**2
    dE_dA = -(centers_sq + 0.5 * inv_A)

    # === Optimized contractions ===
    # Contract over m first: (m,n,n) @ (m,n,n,D) -> (n,n,D), then sum rows/cols
    # This is faster than direct (m,n,n,D) -> (n,D) einsum

    # dE_dA contractions
    dA_ij = jnp.einsum("mij,mijd->ijd", weighted, dE_dA)  # (n, n, D)
    w_dA_row = jnp.sum(dA_ij, axis=1)  # sum over j -> (n, D)
    w_dA_col = jnp.sum(dA_ij, axis=0)  # sum over i -> (n, D)

    # centers contractions
    c_ij = jnp.einsum("mij,mijd->ijd", weighted, centers)  # (n, n, D)
    w_c_row = jnp.sum(c_ij, axis=1)  # (n, D)
    w_c_col = jnp.sum(c_ij, axis=0)  # (n, D)

    # Scalar weights per basis function
    w_ij = jnp.sum(weighted, axis=0)  # (n, n) sum over m
    w_row = jnp.sum(w_ij, axis=1)  # (n,) sum over j
    w_col = jnp.sum(w_ij, axis=0)  # (n,) sum over i

    # === Alpha gradients (symmetric structure) ===
    grad_alpha_i_conj = w_dA_row + 2 * mu * w_c_row - mu**2 * w_row[:, None]
    grad_alpha_j = w_dA_col + 2 * mu * w_c_col - mu**2 * w_col[:, None]

    grad_alpha_combined = grad_alpha_i_conj + grad_alpha_j
    grad_b_antisym = 1j * (grad_alpha_j - grad_alpha_i_conj)

    # === Mu gradients ===
    grad_mu_i = 2 * alpha.conj() * (w_c_row - mu * w_row[:, None]) + 1j * p * w_row[:, None]
    grad_mu_j = 2 * alpha * (w_c_col - mu * w_col[:, None]) - 1j * p * w_col[:, None]
    grad_mu = grad_mu_i + grad_mu_j

    # === P gradients ===
    grad_p_i = 1j * (mu * w_row[:, None] - w_c_row)
    grad_p_j = 1j * (w_c_col - mu * w_col[:, None])
    grad_p = grad_p_i + grad_p_j

    return grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_row + w_col


def assemble_param_gradients_fused(
    grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_sum, a
):
    """
    Assemble final parameter gradients from fused component gradients.

    Args:
        grad_alpha_combined: (n, D) sum of alpha gradients
        grad_b_antisym: (n, D) antisymmetric combination for b gradient
        grad_mu, grad_p: (n, D) combined gradients
        w_sum: (n,) row_sum + col_sum of weighted contributions
        a: (n, D) width parameter

    Returns:
        grad_params: (n, D, 4) real gradients [grad_a, grad_b, grad_mu, grad_p]
    """
    # grad_a = 2*a*(grad_alpha_i_conj + grad_alpha_j) + normalization
    grad_a = 2 * a * grad_alpha_combined + 0.5 / a * w_sum[:, None]

    # grad_b = -1j*grad_alpha_i_conj + 1j*grad_alpha_j (already computed)
    grad_b = grad_b_antisym

    grad_params = jnp.stack([grad_a, grad_b, grad_mu, grad_p], axis=2)
    return jnp.real(grad_params)


# =============================================================================
# RAW COMPUTATION (for naive autodiff reference)
# =============================================================================


def _calculate_Gaussian_expectation_values_raw(
    params, gaussian_exponential_params, linear_params
):
    """
    Raw computation of <g_i|V|g_j> for Gaussian potentials.

    This is the core computation without custom_vjp decoration.
    Used for naive autodiff and internally by other versions.
    """
    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)
    alphaK = aK**2 + 1j * bK  # (m, D)

    allA, allB, allC = combine_with_potential(gauss_A, gauss_B, gauss_C, alphaK, muK, pK)
    contribs, _ = compute_contributions(allA, allB, allC, exp_substr_terms, linear_params)

    return jnp.sum(contribs, axis=0)  # (n, n)


# =============================================================================
# FAST VERSION: Fully vectorized with custom VJP
# =============================================================================


@partial(jax.custom_vjp, nondiff_argnums=(1, 2))
def calculate_Gaussian_expectation_values_fast(
    params, gaussian_exponential_params, linear_params
):
    """
    Fast vectorized computation of <g_i|V|g_j> for Gaussian potentials.

    Vectorizes over all m potential terms for better GPU utilization.
    Uses custom VJP for efficient backward pass.
    """
    return _calculate_Gaussian_expectation_values_raw(
        params, gaussian_exponential_params, linear_params
    )


def _fast_fwd(params, gaussian_exponential_params, linear_params):
    """Forward pass - cache basis terms to avoid recomputation in backward."""
    a, b, mu, p = get_individual_params(params)
    alpha = a**2 + 1j * b

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)
    alphaK = aK**2 + 1j * bK

    allA, allB, allC = combine_with_potential(gauss_A, gauss_B, gauss_C, alphaK, muK, pK)
    contribs, centers = compute_contributions(
        allA, allB, allC, exp_substr_terms, linear_params
    )
    result = jnp.sum(contribs, axis=0)

    # Cache everything needed for backward (avoid recomputation)
    residuals = (a, alpha, mu, p, allA, centers, contribs)
    return result, residuals


def _fast_bwd(gaussian_exponential_params, linear_params, residuals, g):
    """Vectorized backward pass using cached residuals."""
    a, alpha, mu, p, allA, centers, contribs = residuals

    # Weight by output gradient
    weighted = g[None, :, :] * contribs  # (m, n, n)

    # Compute parameter gradients
    grads = compute_param_gradients(weighted, centers, allA, alpha, mu, p)
    grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_sum = grads

    # Assemble final gradients
    grad_params = assemble_param_gradients_fused(
        grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_sum, a
    )
    return (grad_params,)


calculate_Gaussian_expectation_values_fast.defvjp(_fast_fwd, _fast_bwd)


# =============================================================================
# BATCHED VERSION: Configurable memory usage via batching
# =============================================================================


def compute_batch_size(m, n, D, memory_limit_bytes=1_000_000_000):
    """
    Compute batch size based on measured memory scaling.

    Empirically, memory scales as ~1.7 * batch * n * n * D * 16 bytes.
    """
    bytes_per_batch_item = int(1.7 * n * n * D * 16)
    batch_size = max(1, int(memory_limit_bytes / bytes_per_batch_item))
    return min(m, batch_size)


@partial(jax.custom_vjp, nondiff_argnums=(1, 2, 3))
def calculate_Gaussian_expectation_values_batched(
    params, gaussian_exponential_params, linear_params, memory_limit_mb=1000
):
    """
    Batched VJP version that bounds memory usage.

    Args:
        params: (n, D, 4) basis parameters [a, b, mu, p]
        gaussian_exponential_params: (m, D, 4) potential Gaussian parameters
        linear_params: (m,) linear coefficients
        memory_limit_mb: Maximum memory to use in backward pass (default 1GB)

    Returns:
        (n, n) matrix of expectation values
    """
    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)
    alphaK = aK**2 + 1j * bK

    allA, allB, allC = combine_with_potential(gauss_A, gauss_B, gauss_C, alphaK, muK, pK)
    contribs, _ = compute_contributions(allA, allB, allC, exp_substr_terms, linear_params)

    return jnp.sum(contribs, axis=0)


def _batched_fwd(params, gaussian_exponential_params, linear_params, memory_limit_mb):
    """Forward pass storing only params."""
    result = calculate_Gaussian_expectation_values_batched(
        params, gaussian_exponential_params, linear_params, memory_limit_mb
    )
    return result, params


def _batched_bwd(gaussian_exponential_params, linear_params, memory_limit_mb, res, g):
    """
    Batched backward pass that processes potential terms in chunks.

    Memory is bounded by processing only batch_size potential terms at once.
    """
    params = res
    n, D, _ = params.shape
    m = gaussian_exponential_params.shape[0]
    memory_limit_bytes = memory_limit_mb * 1_000_000
    batch_size = compute_batch_size(m, n, D, memory_limit_bytes)

    # Pre-compute basis quantities
    a, b, mu, p = get_individual_params(params)
    alpha = a**2 + 1j * b

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

    # Pad potential params to multiple of batch_size
    num_batches = (m + batch_size - 1) // batch_size
    padded_m = num_batches * batch_size
    pad_amount = padded_m - m

    if pad_amount > 0:
        last_row = gaussian_exponential_params[-1:, :, :]
        pad_exp = jnp.tile(last_row, (pad_amount, 1, 1))
        gaussian_exponential_params_padded = jnp.concatenate(
            [gaussian_exponential_params, pad_exp], axis=0
        )
        pad_lin = jnp.zeros((pad_amount,), dtype=linear_params.dtype)
        linear_params_padded = jnp.concatenate([linear_params, pad_lin], axis=0)
    else:
        gaussian_exponential_params_padded = gaussian_exponential_params
        linear_params_padded = linear_params

    # Reshape into batches
    exp_batched = gaussian_exponential_params_padded.reshape(num_batches, batch_size, D, 4)
    lin_batched = linear_params_padded.reshape(num_batches, batch_size)

    def process_batch(carry, inputs):
        """Process one batch of potential terms."""
        exp_batch, lin_batch = inputs

        aK_batch, bK_batch, muK_batch, pK_batch = exp_batch.transpose(2, 0, 1)
        alphaK_batch = aK_batch**2 + 1j * bK_batch

        allA, allB, allC = combine_with_potential(
            gauss_A, gauss_B, gauss_C, alphaK_batch, muK_batch, pK_batch
        )
        contribs, centers = compute_contributions(
            allA, allB, allC, exp_substr_terms, lin_batch
        )

        weighted = g[None, :, :] * contribs

        grads = compute_param_gradients(weighted, centers, allA, alpha, mu, p)
        # grads = (grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_sum)

        new_carry = jax.tree.map(lambda c, b: c + b, carry, grads)
        return new_carry, None

    # Initialize accumulators
    zero_nD = jnp.zeros((n, D), dtype=jnp.complex128)
    zero_n = jnp.zeros((n,), dtype=jnp.complex128)
    init_carry = (zero_nD, zero_nD, zero_nD, zero_nD, zero_n)

    # Run scan over batches
    final_carry, _ = jax.lax.scan(process_batch, init_carry, (exp_batched, lin_batched))
    grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_sum = final_carry

    # Assemble final gradients
    grad_params = assemble_param_gradients_fused(
        grad_alpha_combined, grad_b_antisym, grad_mu, grad_p, w_sum, a
    )
    return (grad_params,)


calculate_Gaussian_expectation_values_batched.defvjp(_batched_fwd, _batched_bwd)
