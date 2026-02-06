import jax.numpy as jnp


def broadcast_params(params):
    n, D, _ = params.shape
    a = params[:, :, 0]  # (n,D)
    b = params[:, :, 1]
    mu = params[:, :, 2]
    p = params[:, :, 3]
    a_m = jnp.broadcast_to(a[:, None, :], (n, n, D))
    a_n = jnp.broadcast_to(a[None, :, :], (n, n, D))
    b_m = jnp.broadcast_to(b[:, None, :], (n, n, D))
    b_n = jnp.broadcast_to(b[None, :, :], (n, n, D))
    mu_m = jnp.broadcast_to(mu[:, None, :], (n, n, D))
    mu_n = jnp.broadcast_to(mu[None, :, :], (n, n, D))
    p_m = jnp.broadcast_to(p[:, None, :], (n, n, D))
    p_n = jnp.broadcast_to(p[None, :, :], (n, n, D))
    return a_m, a_n, b_m, b_n, mu_m, mu_n, p_m, p_n


def get_terms(params):
    a_m, a_n, b_m, b_n, mu_m, mu_n, p_m, p_n = broadcast_params(params)
    ket_terms = [
        [(-4 * a_m**2 * mu_m**2 + 1) / (2 * a_m), 4 * a_m * mu_m, -2 * a_m],
        [-1j * mu_m**2, 2 * 1j * mu_m, -1j],
        [-2 * a_m**2 * mu_m - 2 * 1j * b_m * mu_m - 1j * p_m, 2 * a_m**2 + 2 * 1j * b_m],
        [-1j * mu_m, 1j],
    ]
    bra_terms = [
        [(-4 * a_n**2 * mu_n**2 + 1) / (2 * a_n), 4 * a_n * mu_n, -2 * a_n],
        [1j * mu_n**2, -2 * 1j * mu_n, 1j],
        [-2 * a_n**2 * mu_n + 2 * 1j * b_n * mu_n + 1j * p_n, 2 * a_n**2 - 2 * 1j * b_n],
        [1j * mu_n, -1j],
    ]
    n, D, _ = params.shape
    max_deg = 3

    zero = jnp.zeros_like(a_m)

    def stack_param_terms(term_list):
        # returns (n, n, D, 4, max_deg)
        coeffs_per_param = []
        for d in range(4):
            pieces = []
            for x in range(max_deg):
                if x < len(term_list[d]):
                    term = term_list[d][x]
                    # broadcast scalars / lower-rank arrays to (n,n,D)
                    term = jnp.broadcast_to(term, zero.shape)
                else:
                    term = zero
                pieces.append(term)
            # stack over degree axis -> (n,n,D,max_deg)
            coeffs_per_param.append(jnp.stack(pieces, axis=-1))
        # stack over parameter axis -> (n,n,D,4,max_deg)
        return jnp.stack(coeffs_per_param, axis=3)

    ket_coeffs = stack_param_terms(ket_terms)
    bra_coeffs = stack_param_terms(bra_terms)
    return ket_coeffs, bra_coeffs


def get_non_diagonal_terms(params, S, moments):
    n, D, _ = params.shape
    M = jnp.zeros((n, n, D, D, 4, 4), dtype=jnp.complex128)
    ket_coeffs, bra_coeffs = get_terms(params)

    # Calculate the expectations for ket terms.
    # The reason we can treat them separately is that they factorize over dimensions for uncorrelated Gaussians.
    m3 = moments[..., :3]  # (n,n,D,3)

    # elementwise multiply coeffs with the correct moment power:
    # ket_expectations[i,j,d,s,x] = ket_coeffs[i,j,d,s,x] * moments[i,j,d,x]
    ket_expectations = jnp.einsum("ijdsx,ijdx->ijdsx", ket_coeffs, m3)
    bra_expectations = jnp.einsum("ijdsx,ijdx->ijdsx", bra_coeffs, m3)
    M = jnp.einsum("ijdsx,ijpty,ij->ijdpst", ket_expectations, bra_expectations, S)
    return M


def build_spd_hessian(params, S_linear_contraction, moments):
    """
    Build a fully SPD Hessian approximation purely from
    ket/bra polynomials + moments (no get_4_times_4_coefficients_*).

    params: (n,D,4)
    S_linear_contraction: (n,n) with entries c_i^* S_ij c_j
    moments: (n,n,D,5)  -- normalized moments up to x^4
    """
    n, D, _ = params.shape
    ket_coeffs, bra_coeffs = get_terms(params)  # (n,n,D,4,3) each

    # moments for degree 0..2 (for separate expectations)
    m3 = moments[..., :3]  # (n,n,D,3)
    # moments for degree 0..4 (for same-dim convolution)
    m5 = moments[..., :5]  # (n,n,D,5)

    # 1) Cross-dimension contribution (d != p)
    # E_ket[i,j,d,s] = sum_x ket_coeffs[i,j,d,s,x] * moments[i,j,d,x]
    E_ket = jnp.einsum("ijdsx,ijdx->ijds", ket_coeffs, m3)  # (n,n,D,4)
    E_bra = jnp.einsum("ijdtx,ijdx->ijdt", bra_coeffs, m3)  # (n,n,D,4)

    # Cross[d,p] = S_ij * E_ket[d,s] * E_bra[p,t]
    full_hessian = jnp.einsum(
        "ijds,ijpt,ij->ijdpst", E_ket, E_bra, S_linear_contraction
    )  # (n,n,D,D,4,4)

    # 2) Same-dimension contribution (d == p) with polynomial product
    # K[i,j,d,s,t] = sum_{k,l=0..2} A_{k}^* B_{l} * m5[k+l]
    K = jnp.zeros((n, n, D, 4, 4), dtype=jnp.complex128)

    for k in range(3):
        for l in range(3):
            mkl = m5[..., k + l]  # (n,n,D)
            K += jnp.einsum(
                "ijd s,ijd t,ijd->ijdst",
                ket_coeffs[..., k],  # (n,n,D,4)
                bra_coeffs[..., l],  # (n,n,D,4)
                mkl,  # (n,n,D)
            )

    # Multiply by S_linear_contraction for the wavefunction coefficients
    same_dim = jnp.einsum("ij,ijdst->ijdst", S_linear_contraction, K)  # (n,n,D,4,4)

    # 3) Insert same-dimension blocks into the full tensor
    idx = jnp.arange(D)
    full_hessian = full_hessian.at[:, :, idx, idx, :, :].set(same_dim)

    # You probably want a real symmetric approximation here
    return full_hessian.real


def calculate_Hessian_matrix(params, S, moments, linear_coefficients):
    """
    Clean SPD Hessian from ket/bra + moments only.
    """
    assert moments.shape[3] >= 5, "Need moments up to x^4 for SPD Hessian calculation."
    # Contract S with linear coefficients: c_i^* S_ij c_j
    S_linear_contraction = jnp.einsum(
        "i,ij,j->ij", linear_coefficients.conj(), S, linear_coefficients
    )  # (n,n)

    H = build_spd_hessian(params, S_linear_contraction, moments)

    return H
