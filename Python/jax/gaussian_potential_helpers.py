import jax
import jax.numpy as jnp
from helpers import get_individual_params

dtype_complex = jnp.complex128
dtype_real = jnp.float64


def get_gaussian_pair_terms(params):
    a, b, mu, p = get_individual_params(params)
    alpha = a**2 + 1j * b  # (n,D)
    ai, aj = alpha.conj()[:, None, :], alpha[None, :, :]
    gauss_A = ai + aj  # (n,n,D)
    ai_mui = ai * (mu[:, None, :])
    ai_mui2 = ai * (mu[:, None, :] ** 2)
    aj_muj = aj * (mu[None, :, :])
    aj_muj2 = aj * (mu[None, :, :] ** 2)
    muipi = mu[:, None, :] * p[:, None, :]
    mujpj = mu[None, :, :] * p[None, :, :]
    muipimujpj = 1j * (-mujpj + muipi)
    pi, pj = p[:, None, :], p[None, :, :]
    pipj = pj - pi
    aimui_ajmuj_sum = ai_mui + aj_muj
    aimuisq_ajmujsq = ai_mui2 + aj_muj2
    gauss_B = 2 * aimui_ajmuj_sum + 1j * pipj
    gauss_C = -aimuisq_ajmujsq + muipimujpj
    return gauss_A, gauss_B, gauss_C


def get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C):
    center_diag = gauss_B.diagonal(axis1=0, axis2=1) / (2 * gauss_A.diagonal(axis1=0, axis2=1))
    constant_diag = gauss_C.diagonal(axis1=0, axis2=1) + center_diag**2 * gauss_A.diagonal(
        axis1=0, axis2=1
    )
    exponent_subtraction = constant_diag - 0.5 * jnp.log(gauss_A.diagonal(axis1=0, axis2=1))
    exponent_subtraction = exponent_subtraction.T
    full_exponent_subtraction = (
        exponent_subtraction[:, None, :] + exponent_subtraction[None, :, :]
    )  # Shape (n,n,D)
    return full_exponent_subtraction


# @jit
def calculate_Gaussian_expectation_values(params, gaussian_exponential_params, linear_params):
    n, D, _ = params.shape  # n Gaussians in D dimensions, _ is always 4
    m, _, _ = gaussian_exponential_params.shape  # m Gaussians in D dimensions, _ is always 4
    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)
    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)  # (m, D)
    alphaK = aK**2 + 1j * bK  # (m, D)

    # A replacement of a "for" loop with jax.lax.fori_loop
    def add_expectation_component(k, acc):
        ak = alphaK[k]  # (D,)
        muk = muK[k]  # (D,)
        pk = pK[k]  # (D,)
        allA = gauss_A + ak
        allB = gauss_B + 2 * (ak * muk) + 1j * pk
        allC = gauss_C - ak * muk**2 - 1j * muk * pk
        new_centers = allB / (2 * allA)
        constants = allC + new_centers**2 * allA
        sum_exponents_all = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms
        sum_exponents_ij = jnp.sum(sum_exponents_all, axis=2)  # (n, n)
        contrib = jnp.exp(sum_exponents_ij) * linear_params[k]  # (n, n)
        return acc + contrib

    init_vals = jnp.zeros((n, n), dtype=dtype_complex)  # The expectation values <g_i|V|g_j>

    expectation_values = jax.lax.fori_loop(0, m, add_expectation_component, init_vals)
    return expectation_values  # Sum over all potential Gaussians


def calculate_squared_Gaussian_potential(
    gaussian_exponential_params, linear_params, tiny=1e-16
):
    # gaussian_exponential_params shape: (m, D, 4) with order [a, b, mu, p] along the last axis
    m, D, _ = gaussian_exponential_params.shape

    new_Gaussian_exponential_params = jnp.zeros((m * (m + 1) // 2, D, 4), dtype=dtype_real)
    new_linear_params = jnp.zeros((m * (m + 1) // 2,), dtype=dtype_complex)
    log_tiny = jnp.log(tiny)
    count = 0
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
            # prune tiny magnitudes safely
            new_c = jnp.where(
                cand_log_real < log_tiny,
                jnp.array(0.0 + 0.0j, dtype=dtype_complex),
                jnp.exp(cand_log),
            )

            # write params (stack per-dimension)
            new_params = jnp.stack([new_a, new_b, new_mu, new_p], axis=-1).astype(dtype_real)
            new_Gaussian_exponential_params = new_Gaussian_exponential_params.at[count].set(
                new_params
            )

            # double off-diagonals (i<j), keep diagonals once
            new_linear_params = new_linear_params.at[count].set(new_c if i == j else 2 * new_c)

            count += 1

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
