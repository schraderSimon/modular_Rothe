from functools import partial

import jax
import jax.numpy as jnp
from jax.numpy import broadcast_to

from .helpers import get_individual_params

dtype_complex = jnp.complex128
dtype_real = jnp.float64


# =============================================================================
# CORE HELPER FUNCTIONS
# =============================================================================
def copy_over_i(x_nD):
    n, D = x_nD.shape
    return jnp.broadcast_to(x_nD[None, :, :], (n, n, D))


def copy_over_j(x_nD):
    n, D = x_nD.shape
    return jnp.broadcast_to(x_nD[:, None, :], (n, n, D))


def get_gaussian_pair_terms(params):
    """
    Compute pairwise Gaussian terms for basis overlap.

    Args:
        params: (n, D, 4) basis parameters [a, b, mu, p]

    Returns:
        gauss_A, gauss_B, gauss_C: (n, n, D) tensors defining the pairwise Gaussian
    """
    a_nD, b_nD, mu_nD, p_nD = get_individual_params(params)
    alpha_nD = a_nD**2 + 1j * b_nD  # (n, D)
    ai_nnD = copy_over_j(alpha_nD.conj())
    aj_nnD = copy_over_i(alpha_nD)  # (n, n, D)
    gauss_A_nnD = ai_nnD + aj_nnD
    co_j_mu_nnD = copy_over_j(mu_nD)
    co_i_mu_nnD = copy_over_i(mu_nD)
    co_i_p_nnD = copy_over_i(p_nD)
    co_j_p_nnD = copy_over_j(p_nD)
    gauss_B_nnD = 2 * (ai_nnD * co_j_mu_nnD + aj_nnD * co_i_mu_nnD) + 1j * (
        co_i_p_nnD - co_j_p_nnD
    )
    gauss_C_nnD = -(ai_nnD * co_j_mu_nnD**2 + aj_nnD * co_i_mu_nnD**2) + 1j * (
        co_j_mu_nnD * co_j_p_nnD - co_i_mu_nnD * co_i_p_nnD
    )
    return gauss_A_nnD, gauss_B_nnD, gauss_C_nnD


def get_exponent_subtraction_terms(gauss_A_nnD, gauss_B_nnD, gauss_C_nnD):
    """Compute normalization terms from basis overlap."""
    n, _, D = gauss_A_nnD.shape
    idx = jnp.arange(n)
    A_diag_Dn = gauss_A_nnD[idx, idx, :].T  # (D, n)
    B_diag_Dn = gauss_B_nnD[idx, idx, :].T
    C_diag_Dn = gauss_C_nnD[idx, idx, :].T

    center_diag_Dn = B_diag_Dn / (2 * A_diag_Dn)
    constant_diag_Dn = C_diag_Dn + center_diag_Dn**2 * A_diag_Dn
    exponent_subtraction_Dn = (constant_diag_Dn - 0.5 * jnp.log(A_diag_Dn)).T  # (n, D)
    returnval_nnD = copy_over_i(exponent_subtraction_Dn) + copy_over_j(exponent_subtraction_Dn)
    return returnval_nnD


def expand_over_nn(x):
    """
    x: (..., D)
    returns: (..., 1, 1, D)
    """
    return x[..., None, None, :]


def add_batch_dim(x):
    return x[None, ...]


def combine_with_potential(gauss_A_nnD, gauss_B_nnD, gauss_C_nnD, alphaK_BD, muK_BD, pK_BD):
    """
    Combine basis pair terms with potential Gaussian parameters.
    B is the batch size.
    gauss_* : (n, n, D)
    alphaK, muK, pK : (B, D)

    Returns:
        allA, allB, allC : (B, n, n, D)
    """
    aK_B11D = expand_over_nn(alphaK_BD)
    mK_B11D = expand_over_nn(muK_BD)
    pK_B11D = expand_over_nn(pK_BD)
    gauss_A_1nnD = add_batch_dim(gauss_A_nnD)  # (1, n, n, D)
    gauss_B_1nnD = add_batch_dim(gauss_B_nnD)
    gauss_C_1nnD = add_batch_dim(gauss_C_nnD)
    # Broadcasting: (1 n n D) + (B 1 1 D) -> (B n n D)
    allA_BnnD = gauss_A_1nnD + aK_B11D
    allB_BnnD = gauss_B_1nnD + 2 * aK_B11D * mK_B11D + 1j * pK_B11D
    allC_BnnD = gauss_C_1nnD - aK_B11D * mK_B11D**2 - 1j * mK_B11D * pK_B11D

    return allA_BnnD, allB_BnnD, allC_BnnD


def compute_contributions(
    allA_BnnD, allB_BnnD, allC_BnnD, exp_substr_terms_nnD, linear_params_B
):
    """
    Compute weighted Gaussian contributions.

    Args:
        allA, allB, allC: (B, n, n, D) combined Gaussian terms
        exp_substr_terms: (n, n, D) normalization terms
        linear_params: (B,) linear coefficients

    Returns:
        contribs: (B, n, n) weighted contributions per potential term
        centers: (B, n, n, D) Gaussian centers (needed for gradients)
    """
    centers_BnnD = allB_BnnD / (2 * allA_BnnD)
    constants_BnnD = allC_BnnD + centers_BnnD**2 * allA_BnnD
    sum_exponents = (
        constants_BnnD - 0.5 * jnp.log(allA_BnnD) - 0.5 * add_batch_dim(exp_substr_terms_nnD)
    )
    sum_exponents_Bnn = jnp.sum(sum_exponents, axis=-1)  # (B, n, n)

    contribs_Bnn = jnp.exp(sum_exponents_Bnn) * linear_params_B[:, None, None]
    return contribs_Bnn, centers_BnnD


def sum_over_batch(weighted_Bnn, x_BnnD):
    # (B,n,n) · (B,n,n,D) → (n,n,D)
    return jnp.einsum("Bij,Bijd->ijd", weighted_Bnn, x_BnnD)


def get_row_col_contractions(x_nnD):
    # returns (sum over i, sum over j)
    return jnp.sum(x_nnD, axis=0), jnp.sum(x_nnD, axis=1)


def _compute_alpha_and_b_gradients(quants):

    w_dA_row_nD = quants["w_dA_row_nD"]
    w_dA_col_nD = quants["w_dA_col_nD"]
    w_c_row_nD = quants["w_c_row_nD"]
    w_c_col_nD = quants["w_c_col_nD"]
    w_row_nD = quants["w_row_nD"]
    w_col_nD = quants["w_col_nD"]
    mu_nD = quants["mu_nD"]

    grad_alpha_i_conj_nD = w_dA_row_nD + 2 * mu_nD * w_c_row_nD - mu_nD**2 * w_row_nD
    grad_alpha_j_nD = w_dA_col_nD + 2 * mu_nD * w_c_col_nD - mu_nD**2 * w_col_nD

    grad_alpha_combined_nD = grad_alpha_i_conj_nD + grad_alpha_j_nD
    grad_b_antisym_nD = 1j * (grad_alpha_j_nD - grad_alpha_i_conj_nD)

    return grad_alpha_combined_nD, grad_b_antisym_nD


def _compute_mu_gradients(quants):

    w_c_row_nD = quants["w_c_row_nD"]
    w_c_col_nD = quants["w_c_col_nD"]
    w_row_nD = quants["w_row_nD"]
    w_col_nD = quants["w_col_nD"]
    alpha_nD = quants["alpha_nD"]
    mu_nD = quants["mu_nD"]
    p_nD = quants["p_nD"]

    grad_mu_i_nD = 2 * alpha_nD.conj() * (w_c_row_nD - mu_nD * w_row_nD) + 1j * p_nD * w_row_nD
    grad_mu_j_nD = 2 * alpha_nD * (w_c_col_nD - mu_nD * w_col_nD) - 1j * p_nD * w_col_nD

    return grad_mu_i_nD + grad_mu_j_nD


def _compute_p_gradients(quants):

    w_c_row_nD = quants["w_c_row_nD"]
    w_c_col_nD = quants["w_c_col_nD"]
    w_row_nD = quants["w_row_nD"]
    w_col_nD = quants["w_col_nD"]
    mu_nD = quants["mu_nD"]

    grad_p_i_nD = 1j * (mu_nD * w_row_nD - w_c_row_nD)
    grad_p_j_nD = 1j * (w_c_col_nD - mu_nD * w_col_nD)

    return grad_p_i_nD + grad_p_j_nD


def _compute_a_gradients(quants, grad_alpha_combined_nD):
    a_nD = quants["a_nD"]
    n, D = a_nD.shape
    w_sum_n = quants["w_sum_n"]

    return 2 * a_nD * grad_alpha_combined_nD + 0.5 / a_nD * broadcast_to(
        w_sum_n[:, None], (n, D)
    )


def compute_param_gradients(
    weighted_Bnn, centers_BnnD, allA_BnnD, alpha_nD, mu_nD, p_nD, a_nD
):
    """
    Compute gradients w.r.t. basis parameters from weighted contributions.

    Optimization: Contract over B first to get (n,n,D), then sum rows/cols.
    This creates smaller intermediates than direct (B,n,n,D) -> (n,D) contraction.

    Args:
        weighted_Bnn: (B, n, n) weighted output gradients
        centers_BnnD: (B, n, n, D) Gaussian centers
        allA_BnnD: (B, n, n, D) combined A terms
        alpha_nD, mu_nD, p_nD, a_nD: (n, D) basis parameters

    Returns:
        Tuple: (grad_a_nD, grad_b_antisym_nD, grad_mu_nD, grad_p_nD)
    """
    # Core quantities
    n, D = alpha_nD.shape
    inv_A_BnnD = 1.0 / allA_BnnD
    centers_sq_BnnD = centers_BnnD**2
    dE_dA_BnnD = -(centers_sq_BnnD + 0.5 * inv_A_BnnD)

    # === Optimized contractions ===
    # Contract over B first: (B,n,n) @ (B,n,n,D) -> (n,n,D), then sum rows/cols
    # This is faster than direct (B,n,n,D) -> (n,D) einsum

    # dE_dA contractions
    dA_nnD = sum_over_batch(weighted_Bnn, dE_dA_BnnD)  # (n, n, D)
    w_dA_col_nD, w_dA_row_nD = get_row_col_contractions(dA_nnD)  # (n, D) each

    # centers contractions
    c_nnD = sum_over_batch(weighted_Bnn, centers_BnnD)  # (n, n, D)
    w_c_col_nD, w_c_row_nD = get_row_col_contractions(c_nnD)  # (n, D) each

    # Scalar weights per basis function
    w_nn = jnp.sum(weighted_Bnn, axis=0)  # (n, n) sum over B
    w_col_n, w_row_n = get_row_col_contractions(w_nn)  # (n, ) each
    w_col_nD = jnp.broadcast_to(w_col_n[:, None], (n, D))
    w_row_nD = jnp.broadcast_to(w_row_n[:, None], (n, D))
    w_sum_n = w_col_n + w_row_n

    # === Package shared quantities ===
    quants = {
        "w_dA_row_nD": w_dA_row_nD,
        "w_dA_col_nD": w_dA_col_nD,
        "w_c_row_nD": w_c_row_nD,
        "w_c_col_nD": w_c_col_nD,
        "w_row_nD": w_row_nD,
        "w_col_nD": w_col_nD,
        "alpha_nD": alpha_nD,
        "mu_nD": mu_nD,
        "p_nD": p_nD,
        "a_nD": a_nD,
        "w_sum_n": w_sum_n,
    }

    # === Compute gradients using helper functions ===
    grad_alpha_combined_nD, grad_b_antisym_nD = _compute_alpha_and_b_gradients(quants)
    grad_mu_nD = _compute_mu_gradients(quants)
    grad_p_nD = _compute_p_gradients(quants)
    grad_a_nD = _compute_a_gradients(quants, grad_alpha_combined_nD)

    return grad_a_nD, grad_b_antisym_nD, grad_mu_nD, grad_p_nD


def assemble_param_gradients_fused(grad_a, grad_b_antisym, grad_mu, grad_p):

    grad_params_nD4 = jnp.stack([grad_a, grad_b_antisym, grad_mu, grad_p], axis=2)
    return jnp.real(grad_params_nD4)


def compute_batch_size(m, n, D, memory_limit_bytes=1_000_000_000):
    """
    Compute batch size based on measured memory scaling.

    Empirically, memory scales as ~1.7 * batch * n * n * D * 16 bytes.
    """
    bytes_per_batch_item = int(1.7 * n * n * D * 16)
    batch_size = max(1, int(memory_limit_bytes / bytes_per_batch_item))
    return min(m, batch_size)


def _pad_and_batch_potential_params(
    gaussian_exponential_params_mD4, linear_params, batch_size, D
):
    """Pad and reshape potential params into equal-sized batches for scan."""
    m = gaussian_exponential_params_mD4.shape[0]
    num_batches = (m + batch_size - 1) // batch_size
    pad_amount = num_batches * batch_size - m

    if pad_amount > 0:
        pad_exp = jnp.tile(gaussian_exponential_params_mD4[-1:, :, :], (pad_amount, 1, 1))
        gaussian_exponential_params_mD4 = jnp.concatenate(
            [gaussian_exponential_params_mD4, pad_exp], axis=0
        )
        pad_lin = jnp.zeros((pad_amount,), dtype=linear_params.dtype)
        linear_params = jnp.concatenate([linear_params, pad_lin], axis=0)

    exp_batched = gaussian_exponential_params_mD4.reshape(num_batches, batch_size, D, 4)
    lin_batched = linear_params.reshape(num_batches, batch_size)
    return exp_batched, lin_batched


def _unpack_and_combine(
    exp_batch_BD4, lin_batch_B, gauss_A_nnD, gauss_B_nnD, gauss_C_nnD, exp_substr_terms_nnD
):
    """Unpack a batch of potential params and compute contributions.

    Returns:
        contribs_Bnn: (B, n, n)
        centers_BnnD: (B, n, n, D)
        allA_BnnD: (B, n, n, D)
    """
    aK_BD, bK_BD, muK_BD, pK_BD = exp_batch_BD4.transpose(2, 0, 1)
    alphaK_BD = aK_BD**2 + 1j * bK_BD
    allA_BnnD, allB_BnnD, allC_BnnD = combine_with_potential(
        gauss_A_nnD, gauss_B_nnD, gauss_C_nnD, alphaK_BD, muK_BD, pK_BD
    )
    contribs_Bnn, centers_BnnD = compute_contributions(
        allA_BnnD, allB_BnnD, allC_BnnD, exp_substr_terms_nnD, lin_batch_B
    )
    return contribs_Bnn, centers_BnnD, allA_BnnD


@partial(jax.custom_vjp, nondiff_argnums=(1, 2, 3))
def calculate_Gaussian_expectation_values_batched(
    params, gaussian_exponential_params_mD4, linear_params, memory_limit_mb=1000
):
    """
    Batched computation of Gaussian expectation values with bounded memory.

    Both forward and backward passes process potential terms in chunks
    of size determined by memory_limit_mb.

    Args:
        params: (n, D, 4) basis parameters [a, b, mu, p]
        gaussian_exponential_params_mD4: (m, D, 4) potential Gaussian parameters
        linear_params: (m,) linear coefficients
        memory_limit_mb: Maximum memory per batch (default 1GB)

    Returns:
        (n, n) matrix of expectation values
    """
    n, D, _ = params.shape
    m = gaussian_exponential_params_mD4.shape[0]
    memory_limit_bytes = memory_limit_mb * 1_000_000
    batch_size = compute_batch_size(m, n, D, memory_limit_bytes)

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

    exp_batched, lin_batched = _pad_and_batch_potential_params(
        gaussian_exponential_params_mD4, linear_params, batch_size, D
    )

    def process_fwd_batch(carry, inputs):
        exp_batch_BD4, lin_batch_B = inputs
        contribs_Bnn, _, _ = _unpack_and_combine(
            exp_batch_BD4, lin_batch_B, gauss_A, gauss_B, gauss_C, exp_substr_terms
        )
        return carry + jnp.sum(contribs_Bnn, axis=0), None

    init = jnp.zeros((n, n), dtype=dtype_complex)
    result, _ = jax.lax.scan(process_fwd_batch, init, (exp_batched, lin_batched))
    return result


def _batched_fwd(params, gaussian_exponential_params_mD4, linear_params, memory_limit_mb):
    """Forward pass storing only params."""
    result = calculate_Gaussian_expectation_values_batched(
        params, gaussian_exponential_params_mD4, linear_params, memory_limit_mb
    )
    return result, params


def _batched_bwd(gaussian_exponential_params_mD4, linear_params, memory_limit_mb, params, g):
    """Batched backward pass reusing shared batching helpers."""
    n, D, _ = params.shape
    m = gaussian_exponential_params_mD4.shape[0]
    memory_limit_bytes = memory_limit_mb * 1_000_000
    batch_size = compute_batch_size(m, n, D, memory_limit_bytes)

    a, b, mu, p = get_individual_params(params)
    alpha = a**2 + 1j * b

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

    exp_batched, lin_batched = _pad_and_batch_potential_params(
        gaussian_exponential_params_mD4, linear_params, batch_size, D
    )

    def process_bwd_batch(carry, inputs):
        exp_batch_BD4, lin_batch_B = inputs
        contribs_Bnn, centers_BnnD, allA_BnnD = _unpack_and_combine(
            exp_batch_BD4, lin_batch_B, gauss_A, gauss_B, gauss_C, exp_substr_terms
        )
        weighted_Bnn = g[None, :, :] * contribs_Bnn
        grads = compute_param_gradients(weighted_Bnn, centers_BnnD, allA_BnnD, alpha, mu, p, a)
        return jax.tree.map(lambda c, b: c + b, carry, grads), None

    zero_nD = jnp.zeros((n, D), dtype=jnp.complex128)
    init_carry = (zero_nD, zero_nD, zero_nD, zero_nD)

    final_carry, _ = jax.lax.scan(process_bwd_batch, init_carry, (exp_batched, lin_batched))
    grad_a, grad_b_antisym, grad_mu, grad_p = final_carry

    grad_params = assemble_param_gradients_fused(grad_a, grad_b_antisym, grad_mu, grad_p)
    return (grad_params,)


calculate_Gaussian_expectation_values_batched.defvjp(_batched_fwd, _batched_bwd)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def calculate_squared_Gaussian_potential(
    gaussian_exponential_params_mD4, linear_params, tiny=1e-16
):
    # gaussian_exponential_params_mD4 shape: (m, D, 4) with order [a, b, mu, p] along the last axis
    m, D, _ = gaussian_exponential_params_mD4.shape

    log_tiny = jnp.log(tiny)

    # Collect valid (non-tiny) terms in lists
    valid_params_list = []
    valid_linear_list = []

    for i in range(m):
        # params for i (each is shape (D,))
        a1, b1, mu1, p1 = gaussian_exponential_params_mD4[i].T
        c1 = linear_params[i]

        a1sq = a1**2
        for j in range(i, m):
            # params for j
            a2, b2, mu2, p2 = gaussian_exponential_params_mD4[j].T
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
