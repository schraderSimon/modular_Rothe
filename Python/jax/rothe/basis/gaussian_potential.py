"""
Batched computation of <g_i|V|g_j> for Gaussian potentials.

Provides a custom VJP so that the backward pass recomputes only what it
needs (centers, allA) from stored basis parameters, keeping peak memory
bounded by ``memory_limit_mb``.
"""

from functools import partial

import jax
import jax.numpy as jnp

dtype_complex = jnp.complex128
dtype_real = jnp.float64


# ---------------------------------------------------------------------------
# Parameter unpacking (inlined from the former helpers.py)
# ---------------------------------------------------------------------------
def get_individual_params(params):
    """Unpack (n, D, 4) params into (a, b, mu, p), each of shape (n, D)."""
    a, b, mu, p = params.transpose(2, 0, 1)
    return a, b, mu, p


# ---------------------------------------------------------------------------
# Core helper functions
# ---------------------------------------------------------------------------
def add_bra_dim(x_kD, b):
    k, D = x_kD.shape
    return jnp.broadcast_to(x_kD[None, :, :], (b, k, D))


def add_ket_dim(x_bD, k):
    b, D = x_bD.shape
    return jnp.broadcast_to(x_bD[:, None, :], (b, k, D))


def get_gaussian_pair_terms(params_ket, params_bra):
    """
    Compute pairwise Gaussian terms for basis overlap.

    Args:
        params: (n, D, 4) basis parameters [a, b, mu, p]

    Returns:
        gauss_A, gauss_B, gauss_C: (n, n, D) tensors defining the pairwise Gaussian
    """
    a_ket_kD, b_ket_kD, mu_ket_kD, p_ket_kD = get_individual_params(params_ket)
    a_bra_bD, b_bra_bD, mu_bra_bD, p_bra_bD = get_individual_params(params_bra)

    alpha_bra_bD = a_bra_bD**2 + 1j * b_bra_bD
    alpha_ket_kD = a_ket_kD**2 + 1j * b_ket_kD

    b = params_bra.shape[0]
    k = params_ket.shape[0]

    alpha_bra_conj_bkD = add_ket_dim(alpha_bra_bD.conj(), k)
    alpha_ket_bkD = add_bra_dim(alpha_ket_kD, b)
    gauss_A_bkD = alpha_bra_conj_bkD + alpha_ket_bkD

    mu_bra_bkD = add_ket_dim(mu_bra_bD, k)
    mu_ket_bkD = add_bra_dim(mu_ket_kD, b)
    p_bra_bkD = add_ket_dim(p_bra_bD, k)
    p_ket_bkD = add_bra_dim(p_ket_kD, b)

    gauss_B_bkD = 2 * (alpha_bra_conj_bkD * mu_bra_bkD + alpha_ket_bkD * mu_ket_bkD) + 1j * (p_ket_bkD - p_bra_bkD)
    gauss_C_bkD = -(alpha_bra_conj_bkD * mu_bra_bkD**2 + alpha_ket_bkD * mu_ket_bkD**2) + 1j * (
        mu_bra_bkD * p_bra_bkD - mu_ket_bkD * p_ket_bkD
    )
    return gauss_A_bkD, gauss_B_bkD, gauss_C_bkD


def get_exponent_subtraction_terms(params_ket, params_bra):
    a_ket = params_ket[:, :, 0]
    a_bra = params_bra[:, :, 0]
    log_normalization_ket_kD = jnp.log(2 * a_ket**2)
    log_normalization_bra_bD = jnp.log(2 * a_bra**2)
    returnval_bkD = -0.5 * (
        add_ket_dim(log_normalization_bra_bD, params_ket.shape[0])
        + add_bra_dim(log_normalization_ket_kD, params_bra.shape[0])
    )
    return returnval_bkD


def expand_over_bk(x):
    """x: (..., D) -> (..., 1, 1, D)"""
    return x[..., None, None, :]


def add_batch_dim(x):
    return x[None, ...]


def combine_with_potential(gauss_A_bkD, gauss_B_bkD, gauss_C_bkD, alphaK_BD, muK_BD, pK_BD):
    """
    Combine basis pair terms with potential Gaussian parameters.

    B is the batch size.
    """
    aK_B11D = expand_over_bk(alphaK_BD)
    mK_B11D = expand_over_bk(muK_BD)
    pK_B11D = expand_over_bk(pK_BD)
    gauss_A_1bkD = add_batch_dim(gauss_A_bkD)
    gauss_B_1bkD = add_batch_dim(gauss_B_bkD)
    gauss_C_1bkD = add_batch_dim(gauss_C_bkD)
    allA_BbkD = gauss_A_1bkD + aK_B11D
    allB_BbkD = gauss_B_1bkD + 2 * aK_B11D * mK_B11D + 1j * pK_B11D
    allC_BbkD = gauss_C_1bkD - aK_B11D * mK_B11D**2 - 1j * mK_B11D * pK_B11D
    return allA_BbkD, allB_BbkD, allC_BbkD


def compute_contributions(allA_BbkD, allB_BbkD, allC_BbkD, exp_substr_terms_bkD, gaussian_linear_params_B):
    """
    Compute weighted Gaussian contributions.

    Returns:
        contribs: (B, n, n) weighted contributions per potential term
        centers: (B, n, n, D) Gaussian centers (needed for gradients)
    """
    centers_BbkD = allB_BbkD / (2 * allA_BbkD)
    constants_BbkD = allC_BbkD + centers_BbkD**2 * allA_BbkD
    sum_exponents = constants_BbkD - 0.5 * jnp.log(allA_BbkD) - 0.5 * add_batch_dim(exp_substr_terms_bkD)
    sum_exponents_Bbk = jnp.sum(sum_exponents, axis=-1)
    contribs_Bbk = jnp.exp(sum_exponents_Bbk) * gaussian_linear_params_B[:, None, None]
    return contribs_Bbk, centers_BbkD


# ---------------------------------------------------------------------------
# Gradient helpers
# ---------------------------------------------------------------------------
def sum_over_batch(weighted_Bbk, x_BbkD):
    return jnp.einsum("Bij,Bijd->ijd", weighted_Bbk, x_BbkD)


def get_row_col_contractions(x_bkD):
    return jnp.sum(x_bkD, axis=0), jnp.sum(x_bkD, axis=1)


# ---------------------------------------------------------------------------
# Batching helpers
# ---------------------------------------------------------------------------
def compute_batch_size(m, b, k, D, memory_limit_bytes=1_000_000_000):
    """Compute batch size based on measured memory scaling."""
    bytes_per_batch_item = int(1.7 * b * k * D * 16)
    batch_size = max(1, int(memory_limit_bytes / bytes_per_batch_item))
    return min(m, batch_size)


def _pad_and_batch_potential_params(gaussian_exponential_params_mD4, gaussian_linear_params, batch_size, D):
    """Pad and reshape potential params into equal-sized batches for scan."""
    m = gaussian_exponential_params_mD4.shape[0]
    num_batches = (m + batch_size - 1) // batch_size
    pad_amount = num_batches * batch_size - m

    if pad_amount > 0:
        pad_exp = jnp.tile(gaussian_exponential_params_mD4[-1:, :, :], (pad_amount, 1, 1))
        gaussian_exponential_params_mD4 = jnp.concatenate([gaussian_exponential_params_mD4, pad_exp], axis=0)
        pad_lin = jnp.zeros((pad_amount,), dtype=gaussian_linear_params.dtype)
        gaussian_linear_params = jnp.concatenate([gaussian_linear_params, pad_lin], axis=0)

    exp_batched = gaussian_exponential_params_mD4.reshape(num_batches, batch_size, D, 4)
    lin_batched = gaussian_linear_params.reshape(num_batches, batch_size)
    return exp_batched, lin_batched


def _unpack_and_combine(exp_batch_BD4, lin_batch_B, gauss_A_bkD, gauss_B_bkD, gauss_C_bkD, exp_substr_terms_bkD):
    """Unpack a batch of potential params and compute contributions."""
    aK_BD, bK_BD, muK_BD, pK_BD = exp_batch_BD4.transpose(2, 0, 1)
    alphaK_BD = aK_BD**2 + 1j * bK_BD
    allA_BbkD, allB_BbkD, allC_BbkD = combine_with_potential(
        gauss_A_bkD, gauss_B_bkD, gauss_C_bkD, alphaK_BD, muK_BD, pK_BD
    )
    contribs_Bbk, centers_BbkD = compute_contributions(
        allA_BbkD, allB_BbkD, allC_BbkD, exp_substr_terms_bkD, lin_batch_B
    )
    return contribs_Bbk, centers_BbkD, allA_BbkD


# ---------------------------------------------------------------------------
# Main function (custom VJP)
# ---------------------------------------------------------------------------
@partial(jax.custom_vjp, nondiff_argnums=(2, 3, 4))
def calculate_Gaussian_expectation_values_batched(
    params_ket, params_bra, gaussian_exponential_params_mD4, gaussian_linear_params, memory_limit_mb=1000
):
    """
    Batched computation of Gaussian expectation values with bounded memory.

    Both forward and backward passes process potential terms in chunks
    of size determined by memory_limit_mb.

    Args:
        params: (n, D, 4) basis parameters [a, b, mu, p]
        gaussian_exponential_params_mD4: (m, D, 4) potential Gaussian parameters
        gaussian_linear_params: (m,) linear coefficients
        memory_limit_mb: Maximum memory per batch (default 1GB)

    Returns:
        (n, n) matrix of expectation values
    """
    k, D, _ = params_ket.shape
    b, D, _ = params_bra.shape
    m = gaussian_exponential_params_mD4.shape[0]
    memory_limit_bytes = memory_limit_mb * 1_000_000
    batch_size = compute_batch_size(m, b, k, D, memory_limit_bytes)

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params_ket, params_bra)
    exp_substr_terms = get_exponent_subtraction_terms(params_ket, params_bra)

    exp_batched, lin_batched = _pad_and_batch_potential_params(
        gaussian_exponential_params_mD4, gaussian_linear_params, batch_size, D
    )

    def process_fwd_batch(carry, inputs):
        exp_batch_BD4, lin_batch_B = inputs
        contribs_Bbk, _, _ = _unpack_and_combine(
            exp_batch_BD4, lin_batch_B, gauss_A, gauss_B, gauss_C, exp_substr_terms
        )
        return carry + jnp.sum(contribs_Bbk, axis=0), None

    init = jnp.zeros((b, k), dtype=dtype_complex)
    result, _ = jax.lax.scan(process_fwd_batch, init, (exp_batched, lin_batched))
    return result


def _batched_fwd(params_ket, params_bra, gaussian_exponential_params_mD4, gaussian_linear_params, memory_limit_mb):
    """Forward pass storing only params."""
    result = calculate_Gaussian_expectation_values_batched(
        params_ket, params_bra, gaussian_exponential_params_mD4, gaussian_linear_params, memory_limit_mb
    )
    return result, (params_ket, params_bra)


def _batched_bwd(gaussian_exponential_params_mD4, gaussian_linear_params, memory_limit_mb, residual, g):
    """Batched backward pass with separate bra/ket gradients.

    Accumulates contracted intermediate quantities (dA, centers, weights)
    over potential-term batches, then assembles the final parameter
    gradients for both ket and bra.
    """
    params_ket, params_bra = residual
    n_ket, D, _ = params_ket.shape
    n_bra = params_bra.shape[0]
    m = gaussian_exponential_params_mD4.shape[0]
    memory_limit_bytes = memory_limit_mb * 1_000_000
    batch_size = compute_batch_size(m, n_bra, n_ket, D, memory_limit_bytes)

    # Unpack ket parameters
    a_ket, b_ket, mu_ket, p_ket = get_individual_params(params_ket)  # each (k, D)
    alpha_ket = a_ket**2 + 1j * b_ket

    # Unpack bra parameters
    a_bra, b_bra, mu_bra, p_bra = get_individual_params(params_bra)  # each (b, D)
    alpha_bra = a_bra**2 + 1j * b_bra

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params_ket, params_bra)
    exp_substr_terms = get_exponent_subtraction_terms(params_ket, params_bra)

    exp_batched, lin_batched = _pad_and_batch_potential_params(
        gaussian_exponential_params_mD4, gaussian_linear_params, batch_size, D
    )

    def process_bwd_batch(carry, inputs):
        exp_batch_BD4, lin_batch_B = inputs
        contribs_Bbk, centers_BbkD, allA_BbkD = _unpack_and_combine(
            exp_batch_BD4, lin_batch_B, gauss_A, gauss_B, gauss_C, exp_substr_terms
        )
        weighted_Bbk = g[None, :, :] * contribs_Bbk  # (B, b, k)

        # Contract over potential batch dimension B -> (b, k, D) / (b, k)
        inv_A_BbkD = 1.0 / allA_BbkD
        dE_dA_BbkD = -(centers_BbkD**2 + 0.5 * inv_A_BbkD)
        dA_bkD = sum_over_batch(weighted_Bbk, dE_dA_BbkD)
        c_bkD = sum_over_batch(weighted_Bbk, centers_BbkD)
        w_bk = jnp.sum(weighted_Bbk, axis=0)

        # Contract to per-ket (sum over bra, axis=0) and per-bra (sum over ket, axis=1)
        batch_result = (
            jnp.sum(dA_bkD, axis=0),  # dA_ket  (k, D)
            jnp.sum(c_bkD, axis=0),  # c_ket   (k, D)
            jnp.sum(w_bk, axis=0),  # w_ket   (k,)
            jnp.sum(dA_bkD, axis=1),  # dA_bra  (b, D)
            jnp.sum(c_bkD, axis=1),  # c_bra   (b, D)
            jnp.sum(w_bk, axis=1),  # w_bra   (b,)
        )
        return jax.tree.map(lambda c, x: c + x, carry, batch_result), None

    init_carry = (
        jnp.zeros((n_ket, D), dtype=jnp.complex128),  # dA_ket
        jnp.zeros((n_ket, D), dtype=jnp.complex128),  # c_ket
        jnp.zeros((n_ket,), dtype=jnp.complex128),  # w_ket
        jnp.zeros((n_bra, D), dtype=jnp.complex128),  # dA_bra
        jnp.zeros((n_bra, D), dtype=jnp.complex128),  # c_bra
        jnp.zeros((n_bra,), dtype=jnp.complex128),  # w_bra
    )

    final_carry, _ = jax.lax.scan(process_bwd_batch, init_carry, (exp_batched, lin_batched))
    dA_ket_kD, c_ket_kD, w_ket_k, dA_bra_bD, c_bra_bD, w_bra_b = final_carry

    # --- Ket gradients ---
    # gauss_A depends on alpha_ket, gauss_B on alpha_ket/mu_ket/p_ket, etc.
    w_ket_kD = jnp.broadcast_to(w_ket_k[:, None], (n_ket, D))
    grad_alpha_ket = dA_ket_kD + 2 * mu_ket * c_ket_kD - mu_ket**2 * w_ket_kD
    grad_mu_ket = 2 * alpha_ket * (c_ket_kD - mu_ket * w_ket_kD) - 1j * p_ket * w_ket_kD
    grad_p_ket = 1j * (c_ket_kD - mu_ket * w_ket_kD)
    # alpha_ket = a_ket^2 + i*b_ket  =>  da: 2a, db: i
    grad_a_ket = 2 * a_ket * grad_alpha_ket + 0.5 / a_ket * w_ket_kD
    grad_b_ket = 1j * grad_alpha_ket

    grad_params_ket = jnp.real(jnp.stack([grad_a_ket, grad_b_ket, grad_mu_ket, grad_p_ket], axis=2))

    # --- Bra gradients ---
    # gauss_A depends on alpha_bra^*, gauss_B on alpha_bra^*/mu_bra/p_bra, etc.
    w_bra_bD = jnp.broadcast_to(w_bra_b[:, None], (n_bra, D))
    grad_alpha_conj_bra = dA_bra_bD + 2 * mu_bra * c_bra_bD - mu_bra**2 * w_bra_bD
    grad_mu_bra = 2 * alpha_bra.conj() * (c_bra_bD - mu_bra * w_bra_bD) + 1j * p_bra * w_bra_bD
    grad_p_bra = -1j * (c_bra_bD - mu_bra * w_bra_bD)
    # alpha_bra^* = a_bra^2 - i*b_bra  =>  da: 2a, db: -i
    grad_a_bra = 2 * a_bra * grad_alpha_conj_bra + 0.5 / a_bra * w_bra_bD
    grad_b_bra = -1j * grad_alpha_conj_bra

    grad_params_bra = jnp.real(jnp.stack([grad_a_bra, grad_b_bra, grad_mu_bra, grad_p_bra], axis=2))

    return (grad_params_ket, grad_params_bra)


calculate_Gaussian_expectation_values_batched.defvjp(_batched_fwd, _batched_bwd)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
def calculate_squared_Gaussian_potential(gaussian_exponential_params_mD4, gaussian_linear_params, tiny=1e-16):
    m, D, _ = gaussian_exponential_params_mD4.shape
    log_tiny = jnp.log(tiny)

    valid_params_list = []
    valid_linear_list = []

    for i in range(m):
        a1, b1, mu1, p1 = gaussian_exponential_params_mD4[i].T
        c1 = gaussian_linear_params[i]
        a1sq = a1**2

        for j in range(i, m):
            a2, b2, mu2, p2 = gaussian_exponential_params_mD4[j].T
            c2 = gaussian_linear_params[j]
            a2sq = a2**2
            A1 = a1sq + 1j * b1
            A2 = a2sq + 1j * b2

            den = a1sq + a2sq
            new_a = jnp.sqrt(den)
            new_b = b1 + b2

            den_pos = den > 0
            mu_gen = (a1sq * mu1 + a2sq * mu2) / jnp.where(den_pos, den, 1.0)
            new_mu = jnp.where(den_pos, mu_gen, 0.5 * (mu1 + mu2))
            new_p = (p1 + p2) + 2 * (b1 * (mu1 - new_mu) + b2 * (mu2 - new_mu))

            per_dim_log = (
                -(A1 * mu1**2 + A2 * mu2**2)
                + (den + 1j * new_b) * new_mu**2
                - 1j * (p1 * mu1 + p2 * mu2)
                + 1j * new_p * new_mu
            )
            log_phase = jnp.sum(per_dim_log)
            cand_log = jnp.log(c1 * c2) + log_phase
            cand_log_real = jnp.real(cand_log)

            if cand_log_real < log_tiny:
                continue

            new_c = jnp.exp(cand_log)
            new_params = jnp.stack([new_a, new_b, new_mu, new_p], axis=-1).astype(dtype_real)
            valid_params_list.append(new_params)
            valid_linear_list.append(new_c if i == j else 2 * new_c)

    if valid_params_list:
        new_Gaussian_exponential_params = jnp.stack(valid_params_list, axis=0)
        new_gaussian_linear_params = jnp.array(valid_linear_list, dtype=dtype_complex)
    else:
        new_Gaussian_exponential_params = jnp.zeros((0, D, 4), dtype=dtype_real)
        new_gaussian_linear_params = jnp.zeros((0,), dtype=dtype_complex)

    return new_Gaussian_exponential_params, new_gaussian_linear_params


def parse_exponential_params(exponential):
    try:
        linear_coefficients, widths, positions = zip(*exponential)
    except ValueError:
        return None, None
    m = len(linear_coefficients)
    D = len(positions[0])
    exponential_params = jnp.zeros((m, D, 4), dtype=dtype_real)
    for k in range(m):
        width = widths[k]
        pos = positions[k]
        a_k = width
        mu_k = jnp.asarray(pos, dtype=dtype_real)
        exponential_params = exponential_params.at[k, :, 0].set(a_k)
        exponential_params = exponential_params.at[k, :, 2].set(mu_k)
    return jnp.asarray(linear_coefficients, dtype=dtype_complex), exponential_params
