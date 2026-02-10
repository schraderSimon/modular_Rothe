# TODOs, separated into small bits:
# DONE # Parse exponential terms into (m, D, 4) arrays
# DONE Calculate Gaussian squared terms
# DONE Include exponential terms into the Hamiltonian implementation
# DONE Include exponential squared terms into the Hamiltonian squared implementation
# DONE Also include cross-terms
# DONE Include time-dependent potentials (I guess this is not too hard)
# DONE Try imaginary time propagation to find ground state instead of energy/variance optimization?
# Include freezing the ground state option
# Include "uniform width" option and make code faster for that case
#
import os

os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true --xla_gpu_autotune_level=4"

import jax

from .helpers import get_individual_params

jax.config.update("jax_enable_x64", True)
jax.config.update("jax_default_matmul_precision", "high")
import jax.numpy as jnp
from jax import jit

from .gaussian_potential_helpers import (
    calculate_Gaussian_expectation_values_batched as calculate_Gaussian_expectation_values,
)
from .gaussian_potential_helpers import (
    calculate_squared_Gaussian_potential,
    get_exponent_subtraction_terms,
    get_gaussian_pair_terms,
    parse_exponential_params,
)

# Memory limit for batched Gaussian potential computation (MB)
GAUSSIAN_POTENTIAL_MEMORY_LIMIT_MB = 100
from .read_string import obtain_polynomial_squared, read_string

logpi = jnp.log(jnp.pi)
log10 = jnp.log(10)

import jax.lax as lax

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

    # phases/prefactors for coefficients
    # kinetic plane-wave phase: exp(-i * sum_d p^2 dt / (2 m_d))
    kin_phase = jnp.exp(1j * jnp.sum(p_nD**2 * dt / (2.0 * masses[None, :]), axis=1))

    gouy_phase = jnp.exp(-0.5j * jnp.sum(jnp.angle(denom_nD), axis=1))
    gamma = kin_phase * gouy_phase  # |gamma| == 1

    new_coeffs = coeffs * gamma
    return new_params, new_coeffs


class ND_potentials:
    """Gaussian-basis building blocks in D dimensions.

    This class precomputes overlap-like quantities and (Hermite) moments for a set
    of (generally non-orthogonal) Gaussian basis functions. Subclasses implement
    a potential model by providing `calculate_H`, `calculate_H2` or `calculate_V`,
    `calculate_V2`.

    Parameter convention (per Gaussian, per dimension):
        - a: width-like parameter (stored as real, used as $a^2$)
        - b: quadratic phase (real)
        - mu: center (real)
        - p: momentum (real)
    Stored in an array of shape (n, D, 4) as (a, b, mu, p).
    """

    def __init__(self, params, K_max=4):
        """
        params: (n, D, 4) -> (a, b, mu, p) per dimension
        K_max:  maximum polynomial order to precompute moments for (static under jit)
        """
        self.params = jnp.asarray(params, dtype=dtype_real)
        self.n = self.params.shape[0]
        self.D = self.params.shape[1]
        self.K_MAX = int(K_max)

        # Unpack per-dimension parameters
        self.a = self.params[:, :, 0]
        self.b = self.params[:, :, 1]
        self.mu = self.params[:, :, 2]
        self.p = self.params[:, :, 3]
        self.S = None
        self.setUpIntermediates()

    def setUpIntermediates(self):
        # shapes: ai,aj,mu_i,mu_j,p_i,p_j are (n,1,D) / (1,n,D)
        alpha_nD = self.a**2 + 1j * self.b  # (n,D)
        ai, aj = alpha_nD.conj()[:, None, :], alpha_nD[None, :, :]  # (n,1,D), (1,n,D)

        self.Gamma = Gamma = ai + aj  # (n,n,D)
        pi = self.p[:, None, :]  # (n,1,D)
        pj = self.p[None, :, :]  # (1,n,D)
        mui = self.mu[:, None, :] - 1j * pi / (2 * ai)
        muj = self.mu[None, :, :] + 1j * pj / (2 * aj)

        self.exponent_contrib = -(pi.conj() ** 2) / (4 * ai) - (pj**2) / (4 * aj)  # (n,n,D)
        self.tilde_k = (ai * aj) / Gamma  # (n,n,D)
        self.tilde_y = mui - muj  # (n,n,D)
        self.tysq = self.tilde_y**2
        self.tksq = self.tilde_k**2

        # Only compute moment-related intermediates if needed (polynomial potentials)
        self.P = (ai * mui + aj * muj) / Gamma  # (n,n,D)
        self.R = 1 / (2 * Gamma)  # (n,n,D)
        self.beta_i = ai / Gamma  # (n,n,D)
        self.beta_j = aj / Gamma  # (n,n,D)
        # Precompute moments up to a fixed order so shapes stay static under jit
        self.moments = self.calculate_all_moments(maximal_order=self.K_MAX)

    def update_parameters(self, params):
        self.params = jnp.asarray(params, dtype=dtype_real)
        self.a = self.params[:, :, 0]
        self.b = self.params[:, :, 1]
        self.mu = self.params[:, :, 2]
        self.p = self.params[:, :, 3]
        self.setUpIntermediates()

    def update_Kmax(self, K_max):
        """Optionally change K_max at runtime (rebuilds cached moments)."""
        self.K_MAX = int(K_max)
        self.moments = self.calculate_all_moments(maximal_order=self.K_MAX)

    def calculate_S(self):
        expo_core = jnp.sum(self.exponent_contrib - self.tysq * self.tilde_k, axis=-1)  # (n,n)
        expo_pref = -0.5 * jnp.sum(jnp.log(self.Gamma), axis=-1)  # (n,n)
        expo = expo_core + expo_pref  # (n,n)

        diag_expo = -0.5 * jnp.sum(jnp.log(2.0 * self.a**2), axis=-1)
        expo_norm = expo - 0.5 * (diag_expo[:, None] + diag_expo[None, :])
        S = jnp.exp(expo_norm)
        # S = jnp.where(jnp.real(expo_norm) < -100, 0.0, jnp.exp(expo_norm))
        self.S = S
        return S

    def calculate_T(self):
        S = self.S if self.S is not None else self.calculate_S()
        T = S * jnp.sum(self.tilde_k - 2 * self.tksq * self.tysq, axis=-1)
        return T

    def calculate_Tsq(self):
        S = self.S if self.S is not None else self.calculate_S()
        F = self.tilde_k - 2 * self.tksq * self.tysq
        G = (
            4 * (self.tksq**2) * (self.tysq**2)
            - 12 * (self.tilde_k**3) * self.tysq
            + 3 * self.tksq
        )
        sumF = jnp.sum(F, axis=-1)
        sumF2 = jnp.sum(F**2, axis=-1)
        sumG = jnp.sum(G, axis=-1)
        Tsq = S * (sumG + sumF**2 - sumF2)
        return Tsq

    def calculate_all_moments(self, maximal_order=None):
        if maximal_order is None:
            maximal_order = self.K_MAX
        K = int(maximal_order)
        P, R = self.P, self.R  # (n,n,D) complex
        M0 = jnp.ones_like(P)
        M1 = P

        if K == 0:
            return M0[..., None].astype(dtype_complex)
        if K == 1:
            return jnp.stack([M0, M1], axis=-1).astype(dtype_complex)

        def body(carry, k):
            M_prev, M_curr = carry
            kf = jnp.asarray(k, dtype=P.real.dtype)  # k = 1..K-1
            M_next = P * M_curr + kf * R * M_prev  # recurrence
            return (M_curr, M_next), M_next

        (_, _), Ms = lax.scan(body, (M0, M1), jnp.arange(1, K))  # Ms shape: (K-1, n, n, D)
        Ms = jnp.moveaxis(Ms, 0, -1)  # -> (n, n, D, K-1)

        moments = jnp.concatenate([M0[..., None], M1[..., None], Ms], axis=-1)  # (n,n,D,K+1)
        return moments.astype(dtype_complex)

    def calculate_polynomial_expectation_value(self, power_string):
        S = self.S if self.S is not None else self.calculate_S()
        power = jnp.asarray(power_string, dtype=jnp.int32)  # (D,)
        moments = self.moments  # (n,n,D,K_MAX+1)

        idx = power[None, None, :, None]
        picked = jnp.take_along_axis(moments, idx, axis=-1)[..., 0]  # (n,n,D)
        result = jnp.prod(picked, axis=-1)  # (n,n)
        return S * result

    def calculate_dipole_moment(self, coeffs):
        """
        Calculate the dipole moment (position expectation value) for all dimensions.

        Args:
            coeffs: Array of shape (n,) with complex linear coefficients

        Returns:
            dipole: Array of shape (D,) with the expectation value <x_d> for each dimension
        """
        S = self.S if self.S is not None else self.calculate_S()
        D = self.D
        dipole = jnp.zeros(D, dtype=dtype_complex)

        for d in range(D):
            # Power string for x_d: [0, 0, ..., 1, ..., 0] with 1 at position d
            power = jnp.zeros(D, dtype=jnp.int32).at[d].set(1)
            x_d_matrix = self.calculate_polynomial_expectation_value(power)
            # <x_d> = c^dagger @ x_d_matrix @ c
            dipole = dipole.at[d].set(jnp.vdot(coeffs, x_d_matrix @ coeffs))

        return jnp.real(dipole)

    def calculate_H(self, t):
        raise NotImplementedError

    def calculate_H2(self, t):
        raise NotImplementedError

    def calculate_SHH2(self, t=0, params=None, splitting_type="none"):
        if params is not None:
            self.update_parameters(params)
        return self._calc_SHH2(t=t, splitting_type=splitting_type)

    # @partial(jax.jit, static_argnames=('splitting_type',), donate_argnums=(0,))
    def _calc_SHH2(
        self, t=0, splitting_type="none"
    ) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
        H: jnp.ndarray
        H2: jnp.ndarray
        S: jnp.ndarray
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
    Polynomial + exponential potential defined by a small tring.

    The `polynomial_string` is parsed by `read_string.read_string` and is expected
    to define the dimensionality and a dict of monomials with coefficients.
    Example (2D):

        dimension 2
        polynomial
        x0x0: 0.5
        x1x1: 0.5
        x0x0x1: 0.111803
        exponential
        1, 0.225345, [0.0, 0.0] # Gaussian with coeff 1, width 0.225345, centered at origin
        0.2, 0.15, [1.0, 0.0] # Gaussian with coeff .2, width 0.15, centered at (1,0)

    Internally, this class also precomputes the squared polynomial (needed for H^2).
    """

    def __init__(self, params, polynomial_string):
        K_max = 1  # Will be increased based on polynomial degree

        self.polynomial, self.exponential = read_string(polynomial_string)
        if list(self.polynomial.keys()) == []:
            self.polynomial = None
        lin_exp, exponent_exp = parse_exponential_params(self.exponential)
        self.linear_exponential_params = lin_exp
        self.exponential_params = exponent_exp
        if self.exponential_params is not None:
            self.exponential_squared_params, self.linear_exponential_squared = (
                calculate_squared_Gaussian_potential(
                    self.exponential_params, self.linear_exponential_params
                )
            )
        else:
            self.exponential_squared_params = None
            self.linear_exponential_squared = None

        if self.polynomial is not None:
            poly_inserted = update_polynomial_values(
                self.polynomial, 0
            )  # Evaluate to zero initially
            polynomial_squared = obtain_polynomial_squared(poly_inserted)
            self._poly_keys = jnp.asarray(list(poly_inserted.keys()), dtype=jnp.int32)

            self._poly_vals = jnp.asarray(list(poly_inserted.values()), dtype=dtype_complex)
            self._poly2_keys = jnp.asarray(list(polynomial_squared.keys()), dtype=jnp.int32)
            self._poly2_vals = jnp.asarray(
                list(polynomial_squared.values()), dtype=dtype_complex
            )
            for key, val in polynomial_squared.items():
                K_max = max(K_max, max(key))
        # Only compute moments if we have polynomial terms
        super().__init__(params, K_max=K_max)
        if self.polynomial is not None:
            self._build_cached_indices()

    def _update_poly_coefficients(self, t):
        # Evaluate all coefficients at time t
        poly_updated = update_polynomial_values(self.polynomial, t)
        self._poly_vals = jnp.asarray(list(poly_updated.values()), dtype=dtype_complex)

        # Recompute squared polynomial (your existing function!)

        poly2 = obtain_polynomial_squared(poly_updated)
        self._poly2_vals = jnp.array(list(poly2.values()))

    def _build_cached_indices(self):
        """Precompute moment-lookup indices for each monomial in the polynomial.

        For a polynomial V = sum_m c_m * prod_d x_d^{n_{m,d}}, we need to look up
        moments M_{n_d}, M_{n_d-1}, M_{n_d-2} from self.moments (shape n,n,D,K+1).

        This builds clipped index arrays (shape D,num_monomials) so that
        jnp.take_along_axis can gather the right moments without Python loops.
        """
        max_moment_order = self.moments.shape[-1] - 1  # K_MAX

        def _make_moment_indices(keys_MD):
            # keys_MD: (num_monomials, D) — power of each dimension per monomial
            powers_DM = jnp.swapaxes(keys_MD, 0, 1)  # (D, num_monomials)
            idx_n_DM = jnp.clip(powers_DM, 0, max_moment_order)
            idx_n1_DM = jnp.clip(powers_DM - 1, 0, max_moment_order)  # for M_{n-1}
            idx_n2_DM = jnp.clip(powers_DM - 2, 0, max_moment_order)  # for M_{n-2}
            return powers_DM, idx_n_DM, idx_n1_DM, idx_n2_DM

        # Indices for the polynomial itself (used in V, VT+TV)
        (
            self._powers_DM,
            self._moment_idx_DM,
            self._moment_idx_n1_DM,
            self._moment_idx_n2_DM,
        ) = _make_moment_indices(self._poly_keys)
        # Indices for the squared polynomial (used in V^2)
        _, self._moment2_idx_DM, _, _ = _make_moment_indices(self._poly2_keys)

    def setUpIntermediates(self):
        super().setUpIntermediates()
        if self.polynomial is not None:
            self._build_cached_indices()

    def calculate_polynomial_V(self, S, t=0) -> jnp.ndarray:
        if self.polynomial is None:
            return jnp.zeros_like(S)
        self._update_poly_coefficients(t)
        M = self._poly_vals.shape[0]
        idx = jnp.broadcast_to(self._moment_idx_DM, self.moments.shape[:-1] + (M,))
        Mn = jnp.take_along_axis(self.moments, idx, axis=-1)  # (n,n,D,M)
        Mprod = jnp.prod(Mn, axis=2)  # (n,n,M)
        return jnp.einsum("ijm,m->ij", S[..., None] * Mprod, self._poly_vals)

    def calculate_gaussian_V(self, S, t=0) -> jnp.ndarray:
        if self.exponential_params is None:
            return jnp.zeros_like(S)
        returnval = calculate_Gaussian_expectation_values(
            self.params,
            self.exponential_params,
            self.linear_exponential_params,
            GAUSSIAN_POTENTIAL_MEMORY_LIMIT_MB,
        )
        return returnval

    def calculate_gaussian_V2(self, S, t=0) -> jnp.ndarray:
        if self.exponential_squared_params is None:
            return jnp.zeros_like(S)
        returnval = calculate_Gaussian_expectation_values(
            self.params,
            self.exponential_squared_params,
            self.linear_exponential_squared,
            GAUSSIAN_POTENTIAL_MEMORY_LIMIT_MB,
        )
        return returnval

    def calculate_V(self, t=0):
        S = self.S if self.S is not None else self.calculate_S()
        pot = self.calculate_polynomial_V(S, t=t)
        pot += self.calculate_gaussian_V(S, t=t)
        return pot

    def calculate_polynomial_V2(self, S, t=0) -> jnp.ndarray:
        if self.polynomial is None:
            return jnp.zeros_like(S)
        self._update_poly_coefficients(t)
        M2 = self._poly2_vals.shape[0]
        idx = jnp.broadcast_to(self._moment2_idx_DM, self.moments.shape[:-1] + (M2,))
        Mn = jnp.take_along_axis(self.moments, idx, axis=-1)
        Mprod = jnp.prod(Mn, axis=2)
        return jnp.einsum("ijm,m->ij", S[..., None] * Mprod, self._poly2_vals)

    def calculate_polynomial_gaussian_cross_terms(self, S, t=0) -> jnp.ndarray:
        # V_poly and V_gauss commute, so the cross term is 2 * <g_i|V_poly V_gauss|g_j>.
        if (
            self.exponential_params is None
            or self.polynomial is None
            or self.linear_exponential_params is None
        ):
            return jnp.zeros_like(S)
        self._update_poly_coefficients(t)
        gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(self.params)
        exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

        lin_params = self.linear_exponential_params
        aK, bK, muK, pK = self.exponential_params.transpose(2, 0, 1)  # (m, D)
        alphaK = aK**2 + 1j * bK  # (m, D)

        dm = self._powers_DM  # (D, num_monomials)
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
            moments = build_moments(Pk, Rk)  # (n,n,D,K_max+1)

            idx = jnp.clip(dm, 0, moments.shape[-1] - 1)
            idx = jnp.broadcast_to(idx, moments.shape[:-1] + (dm.shape[1],))
            Mn = jnp.take_along_axis(moments, idx, axis=-1)
            Mprod = jnp.prod(Mn, axis=2)  # (n,n,M)

            contrib = jnp.einsum("ijm,m->ij", prefactor[..., None] * Mprod, poly_vals)
            return acc + contrib

        m = alphaK.shape[0]
        cross = jax.lax.fori_loop(0, m, add_component, jnp.zeros_like(S, dtype=dtype_complex))
        return 2.0 * cross

    def calculate_V2(self, t=0) -> jnp.ndarray:
        S = self.S if self.S is not None else self.calculate_S()
        polynomial_V2 = self.calculate_polynomial_V2(S, t=t)
        gaussian_V2 = self.calculate_gaussian_V2(S, t=t)
        polynomial_gaussian_cross = self.calculate_polynomial_gaussian_cross_terms(S, t=t)
        return polynomial_V2 + gaussian_V2 + polynomial_gaussian_cross

    def calculate_polynomial_kinetic_cross_terms(self, S_nn, t=0) -> jnp.ndarray:
        """Compute <g_i| (VT + TV) |g_j> for polynomial V, summed over monomials.

        For each monomial x^n = x_1^{n_1} ... x_D^{n_D} with coefficient c_m,
        the VT matrix element decomposes per dimension d as:

            (VT)_d = kinetic_base_d
                   - 2 n_d  k y beta_j  (M_{n_d-1} / M_{n_d})
                   - n_d(n_d-1)/2 beta_j^2  (M_{n_d-2} / M_{n_d})

        where kinetic_base_d = -(2 k^2 y^2 - k) is the pure-kinetic part,
        k = tilde_k, y = tilde_y, beta_j = alpha_j / Gamma, and M_n are the
        Hermite moments of the overlap Gaussian.

        The full VT element is then:
            VT_{ij} = S_{ij} * sum_m c_m * prod_d M_{n_d} * sum_d (VT)_d
        and we return VT + VT^dagger (since TV = VT^dagger for real V).
        """
        if self.polynomial is None:
            return jnp.zeros_like(S_nn, dtype=dtype_complex)
        self._update_poly_coefficients(t)

        # --- Gather overlap intermediates (all shape n,n,D) ---
        k_nnD = self.tilde_k  # alpha_i * alpha_j / Gamma
        y_nnD = self.tilde_y  # tilde_mu_i - tilde_mu_j
        bj_nnD = self.beta_j  # alpha_j / Gamma

        # --- Look up moments M_{n_d}, M_{n_d-1}, M_{n_d-2} for each monomial ---
        M = self._poly_vals.shape[0]  # number of monomials
        broadcast_shape = self.moments.shape[:-1] + (M,)  # (n, n, D, M)

        idx_nnDM = jnp.broadcast_to(self._moment_idx_DM, broadcast_shape)
        idx1_nnDM = jnp.broadcast_to(self._moment_idx_n1_DM, broadcast_shape)
        idx2_nnDM = jnp.broadcast_to(self._moment_idx_n2_DM, broadcast_shape)

        Mn_nnDM = jnp.take_along_axis(self.moments, idx_nnDM, axis=-1)  # M_{n_d}
        Mn1_nnDM = jnp.take_along_axis(self.moments, idx1_nnDM, axis=-1)  # M_{n_d - 1}
        Mn2_nnDM = jnp.take_along_axis(self.moments, idx2_nnDM, axis=-1)  # M_{n_d - 2}

        # Zero out moments for orders below the derivative requirement
        nvec_DM = self._powers_DM  # (D, num_monomials) integer power per dim
        Mn1_nnDM = jnp.where(nvec_DM[None, None, :, :] >= 1, Mn1_nnDM, 0)
        Mn2_nnDM = jnp.where(nvec_DM[None, None, :, :] >= 2, Mn2_nnDM, 0)

        # --- Safe ratios M_{n-1}/M_n and M_{n-2}/M_n ---
        # When M_n == 0 the ratio is irrelevant (multiplied by n_d or n_d*(n_d-1)),
        # so we replace 0 with 1 in the denominator to avoid NaN gradients.
        safe_Mn_nnDM = jnp.where(jnp.abs(Mn_nnDM) > 0, Mn_nnDM, 1.0)
        inv_Mn_nnDM = jnp.where(jnp.abs(Mn_nnDM) > 0, 1.0 / safe_Mn_nnDM, 0.0)
        r1_nnDM = Mn1_nnDM * inv_Mn_nnDM  # M_{n-1} / M_n
        r2_nnDM = Mn2_nnDM * inv_Mn_nnDM  # M_{n-2} / M_n

        # --- Per-dimension VT coefficients (n,n,D,M) ---
        nvec_nnDM = nvec_DM.astype(Mn_nnDM.real.dtype)[None, None, :, :]
        k_nnD1 = k_nnD[..., None]
        y_nnD1 = y_nnD[..., None]
        bj_nnD1 = bj_nnD[..., None]

        kinetic_base_nnD1 = -(2 * k_nnD1**2 * y_nnD1**2 - k_nnD1)
        first_deriv_nnDM = -2 * nvec_nnDM * k_nnD1 * y_nnD1 * bj_nnD1 * r1_nnDM
        second_deriv_nnDM = -0.5 * nvec_nnDM * (nvec_nnDM - 1) * bj_nnD1**2 * r2_nnDM
        vt_coeff_nnDM = kinetic_base_nnD1 + first_deriv_nnDM + second_deriv_nnDM

        # --- Contract: prod over D, sum over D, then weight by polynomial coefficients ---
        monomial_prod_nnM = jnp.prod(Mn_nnDM, axis=2)  # prod_d M_{n_d}
        dim_sum_nnM = jnp.sum(vt_coeff_nnDM, axis=2)  # sum_d (VT)_d
        weighted_nnM = S_nn[..., None] * monomial_prod_nnM * dim_sum_nnM

        VT_nn = jnp.einsum("ijm,m->ij", weighted_nnM, self._poly_vals)
        return VT_nn + VT_nn.conj().T

    def calculate_gaussian_kinetic_cross_terms(self, S, t=0) -> jnp.ndarray:
        """Compute VT for Gaussian potentials (real V), return VT+VT†.

        Vectorized over all m potential terms for GPU efficiency.
        """
        if self.exponential_params is None or self.linear_exponential_params is None:
            return jnp.zeros_like(S, dtype=dtype_complex)

        gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(self.params)
        exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)

        aK, bK, muK, pK = self.exponential_params.transpose(2, 0, 1)
        alphaK = aK**2  # (m, D) real, potentials are real (b=p=0)
        lin_params = self.linear_exponential_params  # (m,)

        alpha_j = (self.a**2 + 1j * self.b)[None, :, :]  # (1,n,D)
        mu_j = self.mu[None, :, :]  # (1,n,D)
        p_j = self.p[None, :, :]  # (1,n,D)

        # Vectorize over m: expand potential params to (m, 1, 1, D)
        ak = alphaK[:, None, None, :]  # (m, 1, 1, D)
        muk = muK[:, None, None, :]  # (m, 1, 1, D)

        # Combine basis pairs with potential: (m, n, n, D)
        allA = gauss_A + ak
        allB = gauss_B + 2 * (ak * muk)
        allC = gauss_C - ak * (muk**2)

        mean = allB / (2 * allA)  # (m, n, n, D)
        var = 1.0 / (2.0 * allA)  # (m, n, n, D)

        constants = allC + mean**2 * allA
        sum_exponents = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms[None, :, :, :]
        prefactor = (
            jnp.exp(jnp.sum(sum_exponents, axis=-1)) * lin_params[:, None, None]
        )  # (m, n, n)

        # Moments relative to potential center (mu_k) and basis center (mu_j)
        E1_k = mean - muk  # (m, n, n, D)
        E2_k = var + E1_k**2

        E_y = mean - mu_j  # (m, n, n, D) - broadcasts (1,n,D) with (m,n,n,D)
        E_y2 = var + E_y**2
        Cov_kj = var + E1_k * E_y

        # VT matrix element terms
        term1 = -0.5 * jnp.sum(4 * (ak**2) * E2_k - 2 * ak, axis=-1)  # (m, n, n)

        term2_core = (-2 * ak) * ((-2 * alpha_j) * Cov_kj + 1j * p_j * E1_k)  # (m, n, n, D)
        term2 = -jnp.sum(term2_core, axis=-1)  # (m, n, n)

        lap_factor = -2 * alpha_j + 4 * (alpha_j**2) * E_y2 - 4j * alpha_j * p_j * E_y - p_j**2
        term3 = -jnp.sum(lap_factor, axis=-1)  # (m, n, n)

        contrib = prefactor * (term1 + term2 + term3)  # (m, n, n)
        VT = jnp.sum(contrib, axis=0)  # (n, n)
        return 0.5 * (VT + VT.conj().T)

    def calculate_TVVT(self, t=0) -> jnp.ndarray:
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

    def calculate_SHH2(self, t=0, params=None, splitting_type="none"):
        S = self.calculate_S()
        if splitting_type == "none":
            H = self.calculate_H(t)
            H2 = self.calculate_H2(t)
        else:
            H = self.calculate_V(t)
            H2 = self.calculate_V2(t)
        return S, H, H2
