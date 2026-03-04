"""
Core physics engine: Gaussian-basis wave function classes.

Defines:
- ``ND_potentials``: base class for overlap, kinetic energy, and moments.
- ``generalPotentialSolver``: subclass for polynomial + Gaussian potentials
  parsed from a string DSL.
- ``propagate_kinetic_analytical``: exact free-particle propagation.
"""

import jax
import jax.lax as lax
import jax.numpy as jnp
from jax import jit
from rothe.basis.gaussian_potential import (
    calculate_Gaussian_expectation_values_batched as calculate_Gaussian_expectation_values,
)
from rothe.basis.gaussian_potential import (
    calculate_squared_Gaussian_potential,
    get_exponent_subtraction_terms,
    get_gaussian_pair_terms,
    get_individual_params,
    parse_exponential_params,
)
from rothe.potential_parser import obtain_polynomial_squared, read_string

# Memory limit for batched Gaussian potential computation (MB)
GAUSSIAN_POTENTIAL_MEMORY_LIMIT_MB = 100

dtype_real = jnp.float64
dtype_complex = jnp.complex128


@jit
def propagate_kinetic_analytical(dt, params, coeffs, masses):
    """
    Exact free-particle (kinetic) propagation for a Gaussian basis.

    Args:
        dt: Time step.
        params: Array of shape (n, D, 4) holding per-Gaussian, per-dimension
            parameters in the order (a, b, mu, p), with $a>0$.
        coeffs: Array of shape (n,) with complex linear coefficients $c_i$.
        masses: Array of shape (D,) with the mass in each dimension.

    Returns:
        (new_params, new_coeffs) after applying $e^{-i T dt}$.
    """
    params = jnp.asarray(params, dtype=dtype_real)
    masses = jnp.asarray(masses, dtype=dtype_real)

    a_nD, b_nD, mu_nD, p_nD = get_individual_params(params)
    alpha_nD = a_nD**2 + 1j * b_nD
    denom_nD = 1.0 + 2j * alpha_nD * dt / masses[None, :]
    alpha_propagated_nD = alpha_nD / denom_nD
    mu_propagated_nD = mu_nD + dt * p_nD / masses[None, :]  # classical drift
    p_propagated_nD = p_nD
    a_propagated_nD = jnp.sqrt(jnp.real(alpha_propagated_nD))
    b_propagated_nD = jnp.imag(alpha_propagated_nD)

    new_params = params.at[:, :, 0].set(a_propagated_nD)
    new_params = new_params.at[:, :, 1].set(b_propagated_nD)
    new_params = new_params.at[:, :, 2].set(mu_propagated_nD)
    new_params = new_params.at[:, :, 3].set(p_propagated_nD)

    kin_phase = jnp.exp(1j * jnp.sum(p_nD**2 * dt / (2.0 * masses[None, :]), axis=1))
    gouy_phase = jnp.exp(-0.5j * jnp.sum(jnp.angle(denom_nD), axis=1))
    gamma = kin_phase * gouy_phase

    new_coeffs = coeffs * gamma
    return new_params, new_coeffs


class ND_potentials:
    """Gaussian-basis building blocks in D dimensions.

    This class precomputes overlap-like quantities and (Hermite) moments for a set
    of (generally non-orthogonal) Gaussian basis functions.  Subclasses implement
    a potential model by providing ``calculate_H``, ``calculate_H2``, etc.

    Parameter convention (per Gaussian, per dimension):
        - a: width-like parameter (stored as real, used as $a^2$)
        - b: quadratic phase (real)
        - mu: center (real)
        - p: momentum (real)

    Stored in an array of shape ``(n, D, 4)`` as ``(a, b, mu, p)``.
    """

    def __init__(self, params_ket, params_bra, K_max=4):
        self.D = params_ket.shape[1]
        self.K_MAX = int(K_max)
        self.set_ket_and_bra_params(params_ket, params_bra)
        self.S = None
        self.setUpIntermediates()

    def set_ket_and_bra_params(self, params_ket, params_bra):
        self.params_ket = jnp.asarray(params_ket, dtype=dtype_real)
        self.params_bra = jnp.asarray(params_bra, dtype=dtype_real)
        self.a_ket = self.params_ket[:, :, 0]
        self.b_ket = self.params_ket[:, :, 1]
        self.mu_ket = self.params_ket[:, :, 2]
        self.p_ket = self.params_ket[:, :, 3]
        self.a_bra = self.params_bra[:, :, 0]
        self.b_bra = self.params_bra[:, :, 1]
        self.mu_bra = self.params_bra[:, :, 2]
        self.p_bra = self.params_bra[:, :, 3]

    def setUpIntermediates(self):
        self.S = None
        # ai, aj = alpha_nD.conj()[:, None, :], alpha_nD[None, :, :]
        ai = (self.a_bra**2 + 1j * self.b_bra).conj()[:, None, :]
        aj = (self.a_ket**2 + 1j * self.b_ket)[None, :, :]
        self.Gamma = Gamma = ai + aj  # Shape: (n_bra, n_ket, D)
        pi = self.p_bra[:, None, :]
        pj = self.p_ket[None, :, :]
        mui = self.mu_bra[:, None, :] - 1j * pi / (2 * ai)
        muj = self.mu_ket[None, :, :] + 1j * pj / (2 * aj)

        self.exponent_contrib = -(pi.conj() ** 2) / (4 * ai) - (pj**2) / (4 * aj)
        self.tilde_k = (ai * aj) / Gamma
        self.tilde_y = mui - muj
        self.tysq = self.tilde_y**2
        self.tksq = self.tilde_k**2

        self.P = (ai * mui + aj * muj) / Gamma
        self.R = 1 / (2 * Gamma)
        self.beta_i = ai / Gamma
        self.beta_j = aj / Gamma
        self.moments = self.calculate_all_moments(maximal_order=self.K_MAX)

    def update_parameters(self, params_ket, params_bra):
        self.set_ket_and_bra_params(params_ket, params_bra)
        self.setUpIntermediates()

    def update_Kmax(self, K_max):
        """Optionally change K_max at runtime (rebuilds cached moments)."""
        self.K_MAX = int(K_max)
        self.moments = self.calculate_all_moments(maximal_order=self.K_MAX)

    def calculate_S(self):
        expo_core = jnp.sum(self.exponent_contrib - self.tysq * self.tilde_k, axis=-1)
        expo_pref = -0.5 * jnp.sum(jnp.log(self.Gamma), axis=-1)
        expo = expo_core + expo_pref

        diag_exp_ket = -0.5 * jnp.sum(jnp.log(2.0 * self.a_ket**2), axis=-1)
        diag_exp_bra = -0.5 * jnp.sum(jnp.log(2.0 * self.a_bra**2), axis=-1)
        expo_norm = expo - 0.5 * (diag_exp_ket[None, :] + diag_exp_bra[:, None])
        S = jnp.exp(expo_norm)
        self.S = S
        return S

    def calculate_T(self):
        S = self.S if self.S is not None else self.calculate_S()
        T = S * jnp.sum(self.tilde_k - 2 * self.tksq * self.tysq, axis=-1)
        return T

    def calculate_Tsq(self):
        S = self.S if self.S is not None else self.calculate_S()
        F = self.tilde_k - 2 * self.tksq * self.tysq
        G = 4 * (self.tksq**2) * (self.tysq**2) - 12 * (self.tilde_k**3) * self.tysq + 3 * self.tksq
        sumF = jnp.sum(F, axis=-1)
        sumF2 = jnp.sum(F**2, axis=-1)
        sumG = jnp.sum(G, axis=-1)
        Tsq = S * (sumG + sumF**2 - sumF2)
        return Tsq

    def calculate_all_moments(self, maximal_order=None, P=None, R=None):
        if maximal_order is None:
            maximal_order = self.K_MAX
        K = int(maximal_order)
        if P is None:
            P = self.P
        if R is None:
            R = self.R
        M0 = jnp.ones_like(P)
        M1 = P

        if K == 0:
            return M0[..., None].astype(dtype_complex)
        if K == 1:
            return jnp.stack([M0, M1], axis=-1).astype(dtype_complex)

        def body(carry, k):
            M_prev, M_curr = carry
            kf = jnp.asarray(k, dtype=P.real.dtype)
            M_next = P * M_curr + kf * R * M_prev
            return (M_curr, M_next), M_next

        (_, _), Ms = lax.scan(body, (M0, M1), jnp.arange(1, K))
        Ms = jnp.moveaxis(Ms, 0, -1)
        moments = jnp.concatenate([M0[..., None], M1[..., None], Ms], axis=-1)
        return moments.astype(dtype_complex)

    def calculate_polynomial_expectation_value(self, power_string):
        S = self.S if self.S is not None else self.calculate_S()
        power = jnp.asarray(power_string, dtype=jnp.int32)
        moments = self.moments
        idx = power[None, None, :, None]
        picked = jnp.take_along_axis(moments, idx, axis=-1)[..., 0]
        result = jnp.prod(picked, axis=-1)
        return S * result

    def calculate_dipole_moment(self, coeffs):
        """
        Calculate the dipole moment (position expectation value) for all dimensions.

        Args:
            coeffs: Array of shape (n,) with complex linear coefficients.

        Returns:
            dipole: Array of shape (D,) with <x_d> for each dimension.
        """
        S = self.S if self.S is not None else self.calculate_S()
        D = self.D
        dipole = jnp.zeros(D, dtype=dtype_complex)

        for d in range(D):
            power = jnp.zeros(D, dtype=jnp.int32).at[d].set(1)
            x_d_matrix = self.calculate_polynomial_expectation_value(power)
            dipole = dipole.at[d].set(jnp.vdot(coeffs, x_d_matrix @ coeffs))

        return jnp.real(dipole)

    def calculate_H(self, t):
        raise NotImplementedError

    def calculate_H2(self, t):
        raise NotImplementedError

    def calculate_SHH2(self, t=0, params_ket=None, params_bra=None, splitting_type="none"):
        if params_bra is None:
            params_bra = params_ket
        if params_ket is not None:
            self.update_parameters(params_ket, params_bra)
        return self._calc_SHH2(t=t, splitting_type=splitting_type)

    def _calc_SHH2(self, t=0, splitting_type="none"):
        S = self.S if self.S is not None else self.calculate_S()
        if splitting_type == "none":
            H = self.calculate_H(t)
            H2 = self.calculate_H2(t)
        else:
            H = self.calculate_V(t)
            H2 = self.calculate_V2(t)
        return S, H, H2

    def calculate_V(self, t):
        raise NotImplementedError

    def calculate_V2(self, t):
        raise NotImplementedError


def update_polynomial_values(polynomial, t):
    if polynomial is None:
        return None
    evaluated_polynomial = {}
    for key, val in polynomial.items():
        evaluated_polynomial[key] = val if not callable(val) else (val(t)).astype(float)
    return evaluated_polynomial


class generalPotentialSolver(ND_potentials):
    """
    Polynomial + exponential potential defined by a string DSL.

    The ``polynomial_string`` is parsed by :func:`rothe.potential_parser.read_string`.
    """

    def __init__(self, params_ket, params_bra, polynomial_string, squared_string=None):
        K_max = 1

        self.polynomial, self.exponential = read_string(polynomial_string)
        if list(self.polynomial.keys()) == []:
            self.polynomial = None
        lin_exp, exponent_exp = parse_exponential_params(self.exponential)
        self.linear_exponential_params = lin_exp
        self.exponential_params = exponent_exp
        if squared_string is not None:
            # Use a pre-fitted approximation of V^2 instead of exact pairwise squaring.
            _, sq_exponential = read_string(squared_string)
            self.linear_exponential_squared, self.exponential_squared_params = parse_exponential_params(sq_exponential)
        elif self.exponential_params is not None:
            self.exponential_squared_params, self.linear_exponential_squared = calculate_squared_Gaussian_potential(
                self.exponential_params, self.linear_exponential_params
            )
        else:
            self.exponential_squared_params = None
            self.linear_exponential_squared = None

        if self.polynomial is not None:
            poly_inserted = update_polynomial_values(self.polynomial, 0)
            polynomial_squared = obtain_polynomial_squared(poly_inserted)
            self._poly_keys = jnp.asarray(list(poly_inserted.keys()), dtype=jnp.int32)
            self._poly_vals = jnp.asarray(list(poly_inserted.values()), dtype=dtype_complex)
            self._poly2_keys = jnp.asarray(list(polynomial_squared.keys()), dtype=jnp.int32)
            self._poly2_vals = jnp.asarray(list(polynomial_squared.values()), dtype=dtype_complex)
            for key, val in polynomial_squared.items():
                K_max = max(K_max, max(key))

        super().__init__(params_ket, params_bra, K_max=K_max)
        if self.polynomial is not None:
            self._build_cached_indices()

    def _update_poly_coefficients(self, t):
        poly_updated = update_polynomial_values(self.polynomial, t)
        self._poly_vals = jnp.asarray(list(poly_updated.values()), dtype=dtype_complex)
        poly2 = obtain_polynomial_squared(poly_updated)
        self._poly2_vals = jnp.array(list(poly2.values()))

    def _build_cached_indices(self):
        """Precompute moment-lookup indices for each monomial in the polynomial."""
        max_moment_order = self.moments.shape[-1] - 1

        def _make_moment_indices(keys_MD):
            powers_DM = jnp.swapaxes(keys_MD, 0, 1)
            idx_n_DM = jnp.clip(powers_DM, 0, max_moment_order)
            idx_n1_DM = jnp.clip(powers_DM - 1, 0, max_moment_order)
            idx_n2_DM = jnp.clip(powers_DM - 2, 0, max_moment_order)
            return powers_DM, idx_n_DM, idx_n1_DM, idx_n2_DM

        (self._powers_DM, self._moment_idx_DM, self._moment_idx_n1_DM, self._moment_idx_n2_DM) = _make_moment_indices(
            self._poly_keys
        )
        _, self._moment2_idx_DM, _, _ = _make_moment_indices(self._poly2_keys)

    def setUpIntermediates(self):
        super().setUpIntermediates()
        if self.polynomial is not None:
            self._build_cached_indices()

    def calculate_polynomial_V(self, S, t=0):
        if self.polynomial is None:
            return jnp.zeros_like(S)
        self._update_poly_coefficients(t)
        M = self._poly_vals.shape[0]
        idx = jnp.broadcast_to(self._moment_idx_DM, self.moments.shape[:-1] + (M,))
        Mn = jnp.take_along_axis(self.moments, idx, axis=-1)
        Mprod = jnp.prod(Mn, axis=2)
        return jnp.einsum("ijm,m->ij", S[..., None] * Mprod, self._poly_vals)

    def calculate_gaussian_V(self, S, t=0):
        if self.exponential_params is None:
            return jnp.zeros_like(S)
        return calculate_Gaussian_expectation_values(
            self.params_ket,
            self.params_bra,
            self.exponential_params,
            self.linear_exponential_params,
            GAUSSIAN_POTENTIAL_MEMORY_LIMIT_MB,
        )

    def calculate_gaussian_V2(self, S, t=0):
        if self.exponential_squared_params is None:
            return jnp.zeros_like(S)
        return calculate_Gaussian_expectation_values(
            self.params_ket,
            self.params_bra,
            self.exponential_squared_params,
            self.linear_exponential_squared,
            GAUSSIAN_POTENTIAL_MEMORY_LIMIT_MB,
        )

    def calculate_V(self, t=0):
        S = self.S if self.S is not None else self.calculate_S()
        pot = self.calculate_polynomial_V(S, t=t)
        pot += self.calculate_gaussian_V(S, t=t)
        return pot

    def calculate_polynomial_V2(self, S, t=0):
        if self.polynomial is None:
            return jnp.zeros_like(S)
        self._update_poly_coefficients(t)
        M2 = self._poly2_vals.shape[0]
        idx = jnp.broadcast_to(self._moment2_idx_DM, self.moments.shape[:-1] + (M2,))
        Mn = jnp.take_along_axis(self.moments, idx, axis=-1)
        Mprod = jnp.prod(Mn, axis=2)
        return jnp.einsum("ijm,m->ij", S[..., None] * Mprod, self._poly2_vals)

    def calculate_polynomial_gaussian_cross_terms(self, S, t=0):
        if self.exponential_params is None or self.polynomial is None or self.linear_exponential_params is None:
            return jnp.zeros_like(S)
        self._update_poly_coefficients(t)
        gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(self.params_ket, self.params_bra)
        exp_substr_terms = get_exponent_subtraction_terms(self.params_ket, self.params_bra)

        lin_params = self.linear_exponential_params
        aK, bK, muK, pK = self.exponential_params.transpose(2, 0, 1)
        alphaK = aK**2 + 1j * bK

        dm = self._powers_DM
        poly_vals = self._poly_vals
        max_order = self.K_MAX

        def build_moments(P, R):
            M0 = jnp.ones_like(P)
            if max_order == 0:
                return M0[..., None].astype(dtype_complex)
            M1 = P
            if max_order == 1:
                return jnp.stack([M0, M1], axis=-1).astype(dtype_complex)

            def body(carry, k):
                M_prev, M_curr = carry
                kf = jnp.asarray(k, dtype=P.real.dtype)
                M_next = P * M_curr + kf * R * M_prev
                return (M_curr, M_next), M_next

            (_, _), Ms = lax.scan(body, (M0, M1), jnp.arange(1, max_order))
            Ms = jnp.moveaxis(Ms, 0, -1)
            moments = jnp.concatenate([M0[..., None], M1[..., None], Ms], axis=-1)
            return moments.astype(dtype_complex)

        def add_component(k, acc):
            ak = alphaK[k]
            muk = muK[k]
            pk = pK[k]
            ck = lin_params[k]

            allA = gauss_A + ak
            allB = gauss_B + 2 * (ak * muk) + 1j * pk
            allC = gauss_C - ak * (muk**2) - 1j * muk * pk

            centers = allB / (2 * allA)
            constants = allC + centers**2 * allA
            sum_exponents = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms
            prefactor = jnp.exp(jnp.sum(sum_exponents, axis=2)) * ck

            Pk = centers
            Rk = 1.0 / (2.0 * allA)
            moments = build_moments(Pk, Rk)

            idx = jnp.clip(dm, 0, moments.shape[-1] - 1)
            idx = jnp.broadcast_to(idx, moments.shape[:-1] + (dm.shape[1],))
            Mn = jnp.take_along_axis(moments, idx, axis=-1)
            Mprod = jnp.prod(Mn, axis=2)

            contrib = jnp.einsum("ijm,m->ij", prefactor[..., None] * Mprod, poly_vals)
            return acc + contrib

        m = alphaK.shape[0]
        cross = jax.lax.fori_loop(0, m, add_component, jnp.zeros_like(S, dtype=dtype_complex))
        return 2.0 * cross

    def calculate_V2(self, t=0):
        S = self.S if self.S is not None else self.calculate_S()
        polynomial_V2 = self.calculate_polynomial_V2(S, t=t)
        gaussian_V2 = self.calculate_gaussian_V2(S, t=t)
        polynomial_gaussian_cross = self.calculate_polynomial_gaussian_cross_terms(S, t=t)
        return polynomial_V2 + gaussian_V2 + polynomial_gaussian_cross

    def _compute_polynomial_VT(self, S_bk, tilde_k, tilde_y, beta_ket, moments):
        """Compute <bra| V*T |ket> for polynomial V, where T acts on ket.

        Args:
            S_bk: Overlap matrix (n_bra, n_ket).
            tilde_k: (n_bra, n_ket, D) Gaussian pair parameter.
            tilde_y: (n_bra, n_ket, D) Gaussian pair parameter.
            beta_ket: alpha_ket / Gamma, shape (n_bra, n_ket, D).
            moments: Hermite moments (n_bra, n_ket, D, K+1).

        Returns:
            VT matrix of shape (n_bra, n_ket).
        """
        M = self._poly_vals.shape[0]
        broadcast_shape = moments.shape[:-1] + (M,)

        idx_DM = jnp.broadcast_to(self._moment_idx_DM, broadcast_shape)
        idx1_DM = jnp.broadcast_to(self._moment_idx_n1_DM, broadcast_shape)
        idx2_DM = jnp.broadcast_to(self._moment_idx_n2_DM, broadcast_shape)

        Mn = jnp.take_along_axis(moments, idx_DM, axis=-1)
        Mn1 = jnp.take_along_axis(moments, idx1_DM, axis=-1)
        Mn2 = jnp.take_along_axis(moments, idx2_DM, axis=-1)

        nvec_DM = self._powers_DM
        Mn1 = jnp.where(nvec_DM[None, None, :, :] >= 1, Mn1, 0)
        Mn2 = jnp.where(nvec_DM[None, None, :, :] >= 2, Mn2, 0)

        safe_Mn = jnp.where(jnp.abs(Mn) > 0, Mn, 1.0)
        inv_Mn = jnp.where(jnp.abs(Mn) > 0, 1.0 / safe_Mn, 0.0)
        r1 = Mn1 * inv_Mn
        r2 = Mn2 * inv_Mn

        nvec = nvec_DM.astype(Mn.real.dtype)[None, None, :, :]
        k1 = tilde_k[..., None]
        y1 = tilde_y[..., None]
        bj1 = beta_ket[..., None]

        kinetic_base = -(2 * k1**2 * y1**2 - k1)
        first_deriv = -2 * nvec * k1 * y1 * bj1 * r1
        second_deriv = -0.5 * nvec * (nvec - 1) * bj1**2 * r2
        vt_coeff = kinetic_base + first_deriv + second_deriv

        monomial_prod = jnp.prod(Mn, axis=2)
        dim_sum = jnp.sum(vt_coeff, axis=2)
        weighted = S_bk[..., None] * monomial_prod * dim_sum

        return jnp.einsum("ijm,m->ij", weighted, self._poly_vals)

    def calculate_polynomial_kinetic_cross_terms(self, S_nn, t=0):
        """Compute <bra_i| (VT + TV) |ket_j> for polynomial V.

        Handles bra != ket by computing VT with the original bra-ket
        pairing, then computing VT with swapped roles to obtain TV
        via TV = conj(VT_swap).T.
        """
        if self.polynomial is None:
            return jnp.zeros_like(S_nn, dtype=dtype_complex)
        self._update_poly_coefficients(t)

        # VT: <bra_i| V*T |ket_j>, T acts on ket
        VT = self._compute_polynomial_VT(S_nn, self.tilde_k, self.tilde_y, self.beta_j, self.moments)

        # TV: <bra_i| T*V |ket_j> = (<ket_j| V*T |bra_i>)*
        # Compute VT with swapped bra<->ket, then conjugate-transpose
        ai_s = (self.a_ket**2 + 1j * self.b_ket).conj()[:, None, :]
        aj_s = (self.a_bra**2 + 1j * self.b_bra)[None, :, :]
        Gamma_s = ai_s + aj_s

        pi_s = self.p_ket[:, None, :]
        pj_s = self.p_bra[None, :, :]
        mui_s = self.mu_ket[:, None, :] - 1j * pi_s / (2 * ai_s)
        muj_s = self.mu_bra[None, :, :] + 1j * pj_s / (2 * aj_s)

        tilde_k_s = (ai_s * aj_s) / Gamma_s
        tilde_y_s = mui_s - muj_s
        beta_j_s = aj_s / Gamma_s  # "ket" in the swapped view is old bra

        P_s = (ai_s * mui_s + aj_s * muj_s) / Gamma_s
        R_s = 1.0 / (2.0 * Gamma_s)
        moments_s = self.calculate_all_moments(P=P_s, R=R_s)

        S_swap = S_nn.conj().T
        VT_swap = self._compute_polynomial_VT(S_swap, tilde_k_s, tilde_y_s, beta_j_s, moments_s)
        TV = VT_swap.conj().T

        return VT + TV

    def calculate_gaussian_kinetic_cross_terms(self, S, t=0):
        """Compute VT+TV for Gaussian potentials (real V)."""
        if self.exponential_params is None or self.linear_exponential_params is None:
            return jnp.zeros_like(S, dtype=dtype_complex)

        gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(self.params_ket, self.params_bra)
        exp_substr_terms = get_exponent_subtraction_terms(self.params_ket, self.params_bra)

        aK, bK, muK, pK = self.exponential_params.transpose(2, 0, 1)
        alphaK = aK**2
        lin_params = self.linear_exponential_params

        alpha_j = (self.a_ket**2 + 1j * self.b_ket)[None, :, :]
        mu_j = self.mu_ket[None, :, :]
        p_j = self.p_ket[None, :, :]

        ak = alphaK[:, None, None, :]
        muk = muK[:, None, None, :]

        allA = gauss_A + ak
        allB = gauss_B + 2 * (ak * muk)
        allC = gauss_C - ak * (muk**2)

        mean = allB / (2 * allA)
        var = 1.0 / (2.0 * allA)

        constants = allC + mean**2 * allA
        sum_exponents = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms[None, :, :, :]
        prefactor = jnp.exp(jnp.sum(sum_exponents, axis=-1)) * lin_params[:, None, None]

        E1_k = mean - muk
        E2_k = var + E1_k**2
        E_y = mean - mu_j
        E_y2 = var + E_y**2
        Cov_kj = var + E1_k * E_y

        term1 = -0.5 * jnp.sum(4 * (ak**2) * E2_k - 2 * ak, axis=-1)
        term2_core = (-2 * ak) * ((-2 * alpha_j) * Cov_kj + 1j * p_j * E1_k)
        term2 = -jnp.sum(term2_core, axis=-1)
        lap_factor = -2 * alpha_j + 4 * (alpha_j**2) * E_y2 - 4j * alpha_j * p_j * E_y - p_j**2
        term3 = -jnp.sum(lap_factor, axis=-1)

        contrib = prefactor * (term1 + term2 + term3)
        VT = jnp.sum(contrib, axis=0)
        return VT

    def calculate_TVVT(self, t=0):
        S = self.S if self.S is not None else self.calculate_S()
        polynomial_kinetic_cross = self.calculate_polynomial_kinetic_cross_terms(S, t=t)
        gaussian_kinetic_cross = self.calculate_gaussian_kinetic_cross_terms(S, t=t)
        return polynomial_kinetic_cross + gaussian_kinetic_cross

    def calculate_H(self, t=0):
        H = self.calculate_T()
        H += self.calculate_V(t)
        return H

    def calculate_H2(self, t=0):
        H2 = self.calculate_Tsq()
        H2 += self.calculate_V2(t)
        H2 += self.calculate_TVVT(t)
        return H2
