"""
Gaussian potential expectation values and their custom VJP implementations.

This module provides several implementations of <g_i|V|g_j> for Gaussian potentials:

1. calculate_Gaussian_expectation_values (scan-based):
   - Memory-efficient: processes one potential term at a time via lax.scan
   - Best for very large m where memory is the bottleneck
   - Custom VJP with analytical gradients

2. calculate_Gaussian_expectation_values_fast:
   - Fully vectorized over all m potential terms
   - Best performance for moderate problem sizes
   - Memory: O(m * n² * D * 16 bytes) per tensor, ~15 tensors in backward pass

3. calculate_Gaussian_expectation_values_batched:
   - Batched version that processes m terms in chunks
   - Configurable memory limit via memory_limit_mb parameter
   - Falls back to fully vectorized when batch_size >= m
   - Speed vs memory trade-off: fewer batches = faster but more memory

Memory estimate for fast/batched:
    mem_gb ≈ 15 * m * n² * D * 16 / 1e9

Also includes:
- calculate_squared_Gaussian_potential: compute V² parameters from V parameters
- parse_exponential_params: parse potential specification into arrays
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
# SCAN-BASED VERSION: Memory-efficient with custom VJP
# =============================================================================


@partial(jax.custom_vjp, nondiff_argnums=(1, 2))
def calculate_Gaussian_expectation_values(params, gaussian_exponential_params, linear_params):
    """
    Compute <g_i|V|g_j> where V is a sum of Gaussian potentials.

    Uses custom VJP for memory-efficient analytical gradients.
    Uses lax.scan to avoid materializing large (m, n, n, D) tensors.

    Args:
        params: Wavefunction parameters, shape (n, D, 4) with [a, b, mu, p]
        gaussian_exponential_params: Potential parameters, shape (m, D, 4)
        linear_params: Potential coefficients, shape (m,)

    Returns:
        Expectation value matrix, shape (n, n)
    """
    n, D, _ = params.shape

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)
    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)
    alphaK = aK**2 + 1j * bK  # (m, D)

    def compute_single_contrib(carry, inputs):
        """Process one potential term at a time."""
        result = carry
        ak, muk, pk, ck = inputs  # Each is (D,) or scalar

        allA = gauss_A + ak  # (n, n, D)
        allB = gauss_B + 2 * (ak * muk) + 1j * pk
        allC = gauss_C - ak * muk**2 - 1j * muk * pk

        new_centers = allB / (2 * allA)
        constants = allC + new_centers**2 * allA
        sum_exponents_all = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms
        sum_exponents_ij = jnp.sum(sum_exponents_all, axis=2)  # (n, n)

        contrib = jnp.exp(sum_exponents_ij) * ck
        return result + contrib, None

    # Stack inputs for scan: (m, D) or (m,)
    scan_inputs = (alphaK, muK, pK, linear_params)
    result, _ = jax.lax.scan(
        compute_single_contrib,
        jnp.zeros((n, n), dtype=dtype_complex),
        scan_inputs,
    )
    return result


def _calculate_Gaussian_expectation_values_fwd(
    params, gaussian_exponential_params, linear_params
):
    """Forward pass: compute result and save minimal residuals for backward pass."""
    n, D, _ = params.shape

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)
    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)
    alphaK = aK**2 + 1j * bK  # (m, D)

    def compute_single_contrib(carry, inputs):
        """Process one potential term, returning contrib for gradient."""
        result = carry
        ak, muk, pk, ck = inputs

        allA = gauss_A + ak
        allB = gauss_B + 2 * (ak * muk) + 1j * pk
        allC = gauss_C - ak * muk**2 - 1j * muk * pk

        new_centers = allB / (2 * allA)
        constants = allC + new_centers**2 * allA
        sum_exponents_all = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms
        sum_exponents_ij = jnp.sum(sum_exponents_all, axis=2)

        contrib = jnp.exp(sum_exponents_ij) * ck
        return result + contrib, contrib  # Return contrib for backward pass

    scan_inputs = (alphaK, muK, pK, linear_params)
    result, all_contribs = jax.lax.scan(
        compute_single_contrib,
        jnp.zeros((n, n), dtype=dtype_complex),
        scan_inputs,
    )
    # all_contribs: (m, n, n)

    # Store minimal residuals - we'll recompute allA and centers in backward
    residuals = (params, all_contribs)
    return result, residuals


def _calculate_Gaussian_expectation_values_bwd(
    gaussian_exponential_params, linear_params, residuals, g
):
    """
    Backward pass: compute analytical gradients w.r.t. params.

    Uses scan to process one potential term at a time, avoiding large tensors.
    """
    params, all_contribs = residuals
    # all_contribs: (m, n, n)

    n, D, _ = params.shape

    a, b, mu, p = get_individual_params(params)
    alpha = a**2 + 1j * b  # (n, D)

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)
    alphaK = aK**2 + 1j * bK  # (m, D)

    # Parameters expanded for broadcasting in gradient computation
    mu_i = mu[:, None, :]  # (n, 1, D)
    p_i = p[:, None, :]
    alpha_i_conj = alpha.conj()[:, None, :]
    mu_j = mu[None, :, :]  # (1, n, D)
    p_j = p[None, :, :]
    alpha_j = alpha[None, :, :]

    def compute_single_grad(carry, inputs):
        """Compute gradient contribution from one potential term."""
        (
            grad_alpha_i_conj,
            grad_alpha_j,
            grad_mu_i,
            grad_mu_j,
            grad_p_i,
            grad_p_j,
            total_weighted,
        ) = carry
        ak, muk, pk, contrib_k = inputs  # ak: (D,), contrib_k: (n, n)

        # Recompute allA and centers for this term
        allA = gauss_A + ak  # (n, n, D)
        allB = gauss_B + 2 * (ak * muk) + 1j * pk
        centers = allB / (2 * allA)  # (n, n, D)

        # Weighted contribution
        weighted = g * contrib_k  # (n, n)
        weighted_3d = weighted[:, :, None]  # (n, n, 1)

        # Gradient components
        inv_A = 1.0 / allA
        dE_dA = -(centers**2 + 0.5 * inv_A)
        dE_dB = centers
        dE_dC = jnp.ones_like(allA)

        # Gradients w.r.t. row parameters (index i)
        dE_dalpha_i_conj = dE_dA + dE_dB * 2 * mu_i + dE_dC * (-(mu_i**2))
        dE_dmu_i = dE_dB * 2 * alpha_i_conj + dE_dC * (-2 * alpha_i_conj * mu_i + 1j * p_i)
        dE_dp_i = dE_dB * (-1j) + dE_dC * (1j * mu_i)

        # Gradients w.r.t. column parameters (index j)
        dE_dalpha_j = dE_dA + dE_dB * 2 * mu_j + dE_dC * (-(mu_j**2))
        dE_dmu_j = dE_dB * 2 * alpha_j + dE_dC * (-2 * alpha_j * mu_j - 1j * p_j)
        dE_dp_j = dE_dB * (1j) + dE_dC * (-1j * mu_j)

        # Accumulate: sum over j for row grads, sum over i for col grads
        grad_alpha_i_conj = grad_alpha_i_conj + jnp.sum(weighted_3d * dE_dalpha_i_conj, axis=1)
        grad_mu_i = grad_mu_i + jnp.sum(weighted_3d * dE_dmu_i, axis=1)
        grad_p_i = grad_p_i + jnp.sum(weighted_3d * dE_dp_i, axis=1)

        grad_alpha_j = grad_alpha_j + jnp.sum(weighted_3d * dE_dalpha_j, axis=0)
        grad_mu_j = grad_mu_j + jnp.sum(weighted_3d * dE_dmu_j, axis=0)
        grad_p_j = grad_p_j + jnp.sum(weighted_3d * dE_dp_j, axis=0)

        total_weighted = total_weighted + weighted

        return (
            grad_alpha_i_conj,
            grad_alpha_j,
            grad_mu_i,
            grad_mu_j,
            grad_p_i,
            grad_p_j,
            total_weighted,
        ), None

    # Initialize accumulators
    init_carry = (
        jnp.zeros((n, D), dtype=dtype_complex),  # grad_alpha_i_conj
        jnp.zeros((n, D), dtype=dtype_complex),  # grad_alpha_j
        jnp.zeros((n, D), dtype=dtype_complex),  # grad_mu_i
        jnp.zeros((n, D), dtype=dtype_complex),  # grad_mu_j
        jnp.zeros((n, D), dtype=dtype_complex),  # grad_p_i
        jnp.zeros((n, D), dtype=dtype_complex),  # grad_p_j
        jnp.zeros((n, n), dtype=dtype_complex),  # total_weighted
    )

    scan_inputs = (alphaK, muK, pK, all_contribs)
    (
        grad_alpha_i_conj,
        grad_alpha_j,
        grad_mu_i,
        grad_mu_j,
        grad_p_i,
        grad_p_j,
        total_weighted,
    ), _ = jax.lax.scan(
        compute_single_grad,
        init_carry,
        scan_inputs,
    )

    # Convert alpha gradients to a and b gradients
    grad_a = grad_alpha_i_conj * 2 * a + grad_alpha_j * 2 * a
    grad_b = grad_alpha_i_conj * (-1j) + grad_alpha_j * 1j
    grad_mu = grad_mu_i + grad_mu_j
    grad_p = grad_p_i + grad_p_j

    # Gradient from normalization
    row_sum = jnp.sum(total_weighted, axis=1)
    col_sum = jnp.sum(total_weighted, axis=0)
    grad_a_norm = 0.5 / a * (row_sum[:, None] + col_sum[:, None])
    grad_a = grad_a + grad_a_norm

    # Stack and take real part
    grad_params = jnp.stack([grad_a, grad_b, grad_mu, grad_p], axis=2)
    grad_params = jnp.real(grad_params)

    return (grad_params,)


calculate_Gaussian_expectation_values.defvjp(
    _calculate_Gaussian_expectation_values_fwd, _calculate_Gaussian_expectation_values_bwd
)


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
# UTILITY FUNCTIONS
# =============================================================================


def calculate_squared_Gaussian_potential(
    gaussian_exponential_params, linear_params, tiny=1e-16
):
    # gaussian_exponential_params shape: (m, D, 4) with order [a, b, mu, p] along the last axis
    m, D, _ = gaussian_exponential_params.shape

    log_tiny = jnp.log(tiny)

    # Collect valid (non-tiny) terms in lists
    valid_params_list = []
    valid_linear_list = []

    for i in range(m):
        # params for i (each is shape (D,))
        a1, b1, mu1, p1 = gaussian_exponential_params[i].T
        c1 = linear_params[i]

        a1sq = a1**2
        for j in range(i, m):
            # params for j
            a2, b2, mu2, p2 = gaussian_exponential_params[j].T
            c2 = linear_params[j]

            a2sq = a2**2
            A1 = a1sq + 1j * b1
            A2 = a2sq + 1j * b2

            den = a1sq + a2sq  # = a_12^2, shape (D,)
            new_a = jnp.sqrt(den)  # real, >=0, per-dimension
            new_b = b1 + b2

            # handle degenerate case den==0 (i.e., a1=a2=0 in that dimension)
            den_pos = den > 0
            mu_gen = (a1sq * mu1 + a2sq * mu2) / jnp.where(den_pos, den, 1.0)
            new_mu = jnp.where(den_pos, mu_gen, 0.5 * (mu1 + mu2))

            # works for both general and degenerate cases
            new_p = (p1 + p2) + 2 * (b1 * (mu1 - new_mu) + b2 * (mu2 - new_mu))

            # amplitude: c12 = (c1*c2) * exp(sum_d [ ... ])
            per_dim_log = (
                -(A1 * mu1**2 + A2 * mu2**2)
                + (den + 1j * new_b) * new_mu**2
                - 1j * (p1 * mu1 + p2 * mu2)
                + 1j * new_p * new_mu
            )  # shape (D,)
            log_phase = jnp.sum(per_dim_log)  # scalar complex

            cand_log = jnp.log(c1 * c2) + (log_phase)  # complex scalar
            cand_log_real = jnp.real(cand_log)

            # Skip tiny magnitudes entirely
            if cand_log_real < log_tiny:
                continue

            new_c = jnp.exp(cand_log)

            # write params (stack per-dimension)
            new_params = jnp.stack([new_a, new_b, new_mu, new_p], axis=-1).astype(dtype_real)
            valid_params_list.append(new_params)

            # double off-diagonals (i<j), keep diagonals once
            valid_linear_list.append(new_c if i == j else 2 * new_c)

    # Stack collected valid terms into arrays
    if valid_params_list:
        new_Gaussian_exponential_params = jnp.stack(valid_params_list, axis=0)
        new_linear_params = jnp.array(valid_linear_list, dtype=dtype_complex)
    else:
        # Return empty arrays with correct shape if no valid terms
        new_Gaussian_exponential_params = jnp.zeros((0, D, 4), dtype=dtype_real)
        new_linear_params = jnp.zeros((0,), dtype=dtype_complex)

    return new_Gaussian_exponential_params, new_linear_params


def parse_exponential_params(exponential):
    try:
        linear_coefficients, widths, positions = zip(*exponential)
    except ValueError as e:
        return None, None
    m = len(linear_coefficients)
    D = len(positions[0])
    exponential_params = jnp.zeros((m, D, 4), dtype=dtype_real)
    for k in range(m):
        width = widths[k]
        pos = positions[k]
        a_k = width
        # a_k = jnp.sqrt(width)  # This is not necessary
        mu_k = jnp.asarray(pos, dtype=dtype_real)
        exponential_params = exponential_params.at[k, :, 0].set(a_k)
        exponential_params = exponential_params.at[k, :, 2].set(mu_k)
    return jnp.asarray(linear_coefficients, dtype=dtype_complex), exponential_params
