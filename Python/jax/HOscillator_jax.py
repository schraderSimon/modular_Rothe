import os
os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true --xla_gpu_autotune_level=4"
from functools import partial
import jax
jax.config.update('jax_enable_x64', True)
jax.config.update('jax_default_matmul_precision', 'high')  # helps on GPU/TPU
import jax.numpy as jnp

from jax.nn import softplus
from jax import jit, grad
from read_string import *
import numpy as np
import jax.scipy as jsp
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.optimize import minimize
logpi = jnp.log(jnp.pi)
log10 = jnp.log(10)
import numpy as np
import jax.lax as lax
from functools import partial
dtype_real=jnp.float64
dtype_complex=jnp.complex128
    
@jit
def propagate_kinetic_analytical(dt, params, lincoeff, masses):
    """
    params: shape (n, d, 4) with (a, b, mu, p) per dim; a>0, b real.
    lincoeff: shape (n,) complex coefficients for each Gaussian.
    masses: shape (d,) masses per dimension.
    normalized:
        True  -> basis Gaussians are normalized after updating (recommended).

    Returns: new_params, new_lincoeff
    """
    params = jnp.asarray(params, dtype=dtype_real)
    masses = jnp.asarray(masses, dtype=dtype_real)

    a = params[:, :, 0]
    b = params[:, :, 1]
    x = params[:, :, 2]
    p = params[:, :, 3]
    alpha = a**2 + 1j*b
    denom = 1.0 + 2j*alpha*dt/masses[None, :]
    alpha_p = alpha / denom
    x_p = x + dt*p/masses[None, :]  # classical drift
    p_p = p
    a_p = jnp.sqrt(jnp.real(alpha_p))
    b_p = jnp.imag(alpha_p)

    new_params = params.at[:, :, 0].set(a_p)
    new_params = new_params.at[:, :, 1].set(b_p)
    new_params = new_params.at[:, :, 2].set(x_p)
    new_params = new_params.at[:, :, 3].set(p_p)

    # phases/prefactors for coefficients
    # kinetic plane-wave phase: exp(-i * sum_d p^2 dt / (2 m_d))
    kin_phase = jnp.exp(1j * jnp.sum(p**2 * dt / (2.0*masses[None, :]), axis=1))

    # We ignore the "normalized"=True - it is always true here. (this is for jax.jit compatibility)
    gouy_phase = jnp.exp(-0.5j * jnp.sum(jnp.angle(denom), axis=1))
    gamma = kin_phase * gouy_phase                             # |gamma| == 1
    
    new_lincoeff = lincoeff * gamma
    return new_params, new_lincoeff
def calculate_squared_Gaussian_potential(gaussian_exponential_params, linear_params, tiny=1e-16):
    # gaussian_exponential_params shape: (m, D, 4) with order [a, b, mu, p] along the last axis
    m, D, _ = gaussian_exponential_params.shape

    new_Gaussian_exponential_params = jnp.zeros((m*(m+1)//2, D, 4), dtype=dtype_real)
    new_linear_params = jnp.zeros((m*(m+1)//2,), dtype=dtype_complex)
    log_tiny=jnp.log(tiny)
    count = 0
    for i in range(m):
        # params for i (each is shape (D,))
        a1 = gaussian_exponential_params[i, :, 0]
        b1 = gaussian_exponential_params[i, :, 1]
        mu1 = gaussian_exponential_params[i, :, 2]
        p1 = gaussian_exponential_params[i, :, 3]
        c1 = linear_params[i]

        a1sq = a1**2
        for j in range(i, m):
            # params for j
            a2 = gaussian_exponential_params[j, :, 0]
            b2 = gaussian_exponential_params[j, :, 1]
            mu2 = gaussian_exponential_params[j, :, 2]
            p2 = gaussian_exponential_params[j, :, 3]
            c2 = linear_params[j]

            a2sq = a2**2
            A1 = a1sq + 1j*b1
            A2 = a2sq + 1j*b2

            den = a1sq + a2sq                      # = a_12^2, shape (D,)
            new_a = jnp.sqrt(den)                  # real, >=0, per-dimension
            new_b = b1 + b2

            # handle degenerate case den==0 (i.e., a1=a2=0 in that dimension)
            den_pos = den > 0
            mu_gen = (a1sq * mu1 + a2sq * mu2) / jnp.where(den_pos, den, 1.0)
            new_mu = jnp.where(den_pos, mu_gen, 0.5*(mu1 + mu2))

            # works for both general and degenerate cases
            new_p = (p1 + p2) + 2*(b1*(mu1 - new_mu) + b2*(mu2 - new_mu))

            # amplitude: c12 = (c1*c2) * exp(sum_d [ ... ])
            per_dim_log = (
                -(A1*mu1**2 + A2*mu2**2)
                + (den + 1j*new_b)*new_mu**2
                - 1j*(p1*mu1 + p2*mu2)
                + 1j*new_p*new_mu
            )  # shape (D,)
            log_phase = jnp.sum(per_dim_log)  # scalar complex

            cand_log = jnp.log(c1 * c2) +(log_phase)  # complex scalar
            cand_log_real = jnp.real(cand_log)
            # prune tiny magnitudes safely
            new_c = jnp.where(cand_log_real < log_tiny, jnp.array(0.0+0.0j, dtype=dtype_complex), jnp.exp(cand_log))

            # write params (stack per-dimension)
            new_params = jnp.stack([new_a, new_b, new_mu, new_p], axis=-1).astype(dtype_real)
            new_Gaussian_exponential_params = new_Gaussian_exponential_params.at[count].set(new_params)

            # double off-diagonals (i<j), keep diagonals once
            new_linear_params = new_linear_params.at[count].set(new_c if i == j else 2*new_c)

            count += 1

    return new_Gaussian_exponential_params, new_linear_params
    
#@jit
def calculate_Gaussian_expectation_values(params,gaussian_exponential_params,linear_params):
    print(params.shape)
    n,D,_ = params.shape # n Gaussians in D dimensions, _ is always 4
    m,_,_ = gaussian_exponential_params.shape # m Gaussians in D dimensions, _ is always 4
    gep=gaussian_exponential_params
    expectation_values = jnp.zeros((n,n),dtype=dtype_complex) # The expectation values <g_i|V|g_j>
    a, b, mu, p = params.transpose(2, 0, 1)
    alpha = a**2 + 1j*b                     # (n,D)
    ai, aj = alpha.conj()[:, None, :], alpha[None, :, :]
    gauss_A = ai + aj                      # (n,n,D)
    ai_mui = ai * (mu[:, None, :])
    ai_mui2= ai * (mu[:, None, :]**2)
    aj_muj = aj * (mu[None, :, :])
    aj_muj2= aj * (mu[None, :, :]**2)
    muipi=mu[:,None,:]*p[:,None,:]
    mujpj=mu[None,:, :]*p[None,:, :]
    muipimujpj=1j*(-mujpj+muipi)
    pi, pj = p[:, None, :], p[None, :, :]
    pipj=pj-pi
    aimui_ajmuj_sum= ai_mui + aj_muj
    aimuisq_ajmujsq= ai_mui2 + aj_muj2
    gauss_B=2*aimui_ajmuj_sum+1j*pipj
    gauss_C= -aimuisq_ajmujsq + muipimujpj
    #Calculate the normalization factor for each Gaussian pair
    center_diag   = gauss_B.diagonal(axis1=0, axis2=1) / (2 * gauss_A.diagonal(axis1=0, axis2=1))
    constant_diag = gauss_C.diagonal(axis1=0, axis2=1) + center_diag**2 * gauss_A.diagonal(axis1=0, axis2=1)
    exponent_subtraction = constant_diag - 0.5 * jnp.log(gauss_A.diagonal(axis1=0, axis2=1))
    center_diag=center_diag.T
    constant_diag=constant_diag.T
    exponent_subtraction=exponent_subtraction.T
    full_exponent_subtraction= exponent_subtraction[:,None,:]+exponent_subtraction[None,:,:] # Shape (n,n,D)
    aK, bK, muK, pK = gep.transpose(2, 0, 1)      # (m, D)
    alphaK = aK**2 + 1j*bK                         # (m, D)

    # replace the entire `for k in range(m): ...` with:
    def _body(k, acc):
        ak  = alphaK[k]        # (D,)
        muk = muK[k]           # (D,)
        pk  = pK[k]            # (D,)
        allA = gauss_A + ak
        allB = gauss_B + 2*(ak*muk) + 1j*pk
        allC = gauss_C - ak*muk**2 - 1j*muk*pk
        new_centers = allB / (2*allA)
        constants   = allC + new_centers**2 * allA
        sum_exponents_all = constants - 0.5*jnp.log(allA) - 0.5*full_exponent_subtraction
        sum_exponents_ij  = jnp.sum(sum_exponents_all, axis=2)   # (n, n)
        contrib = jnp.exp(sum_exponents_ij) * linear_params[k]   # (n, n)
        return acc + contrib

    expectation_values = jax.lax.fori_loop(0, m, _body, expectation_values)
    return expectation_values #Sum over all potential Gaussians
class ND_potentials:
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
        self.a  = self.params[:, :, 0]
        self.b  = self.params[:, :, 1]
        self.mu = self.params[:, :, 2]
        self.p  = self.params[:, :, 3]
        self.S = None
        self.setUpIntermediates()

    def setUpIntermediates(self):
        # shapes: ai,aj,mu_i,mu_j,p_i,p_j are (n,1,D) / (1,n,D)
        alpha = self.a**2 + 1j*self.b                     # (n,D)
        ai, aj = alpha.conj()[:, None, :], alpha[None, :, :]  # (n,1,D), (1,n,D)

        self.Gamma = Gamma = ai + aj                      # (n,n,D)

        pi, pj = self.p[:, None, :], self.p[None, :, :]
        mui = self.mu[:, None, :] - 1j * pi/(2*ai)
        muj = self.mu[None, :, :] + 1j * pj/(2*aj)

        self.exponent_contrib = -(pi.conj()**2)/(4*ai) - (pj**2)/(4*aj)     # (n,n,D)
        self.tilde_k   = (ai * aj) / Gamma                                  # (n,n,D)
        self.tilde_y   = mui - muj                                           # (n,n,D)

        # Optional extras
        self.P       = (ai*mui + aj*muj)/Gamma                               # (n,n,D)
        self.R       = 1/(2*Gamma)                                           # (n,n,D)
        self.beta_i  = ai/Gamma                                              # (n,n,D)
        self.beta_j  = aj/Gamma                                              # (n,n,D)
        self.tysq    = self.tilde_y**2
        self.tksq    = self.tilde_k**2

        # Precompute and cache moments up to a fixed order so shapes stay static under jit
        self.moments = self.calculate_all_moments(maximal_order=self.K_MAX)

    def update_parameters(self, params):
        self.params = jnp.asarray(params, dtype=dtype_real)
        self.a  = self.params[:, :, 0]
        self.b  = self.params[:, :, 1]
        self.mu = self.params[:, :, 2]
        self.p  = self.params[:, :, 3]
        self.setUpIntermediates()

    def update_Kmax(self, K_max):
        """Optionally change K_max at runtime (rebuilds cached moments)."""
        self.K_MAX = int(K_max)
        self.moments = self.calculate_all_moments(maximal_order=self.K_MAX)
    def calculate_S(self):
        D = self.D
        expo_core = jnp.sum(self.exponent_contrib - self.tysq * self.tilde_k, axis=-1)  # (n,n)
        #expo_pref = 0.5 * (D * jnp.log(jnp.pi) - jnp.sum(jnp.log(self.Gamma), axis=-1)) # (n,n) #PI ISN'T NEEDED
        expo_pref = -0.5* jnp.sum(jnp.log(self.Gamma), axis=-1) # (n,n)
        expo = expo_core + expo_pref  # (n,n)

        diag_expo = jnp.real(jnp.diagonal(expo))
        expo_norm = expo - 0.5 * (diag_expo[:, None] + diag_expo[None, :])
        S=jnp.exp(expo_norm)
        #S = jnp.where(jnp.real(expo_norm) < -100, 0.0, jnp.exp(expo_norm))
        self.S = S
        return S

    def calculate_T(self):
        S = self.S if self.S is not None else self.calculate_S()
        T = S * jnp.sum(self.tilde_k - 2*self.tksq*self.tysq, axis=-1)
        return T

    def calculate_Tsq(self):
        S = self.S if self.S is not None else self.calculate_S()
        F = self.tilde_k - 2 * self.tksq * self.tysq
        G = 4 * (self.tksq**2) * (self.tysq**2) - 12 * (self.tilde_k**3) * self.tysq + 3 * self.tksq
        sumF  = jnp.sum(F, axis=-1)
        sumF2 = jnp.sum(F**2, axis=-1)
        sumG  = jnp.sum(G, axis=-1)
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
            kf = jnp.asarray(k, dtype=P.real.dtype)           # k = 1..K-1
            M_next = P * M_curr + kf * R * M_prev            # recurrence
            return (M_curr, M_next), M_next

        (_, _), Ms = lax.scan(body, (M0, M1), jnp.arange(1, K))  # Ms shape: (K-1, n, n, D)
        Ms = jnp.moveaxis(Ms, 0, -1)                             # -> (n, n, D, K-1)

        moments = jnp.concatenate([M0[..., None], M1[..., None], Ms], axis=-1)  # (n,n,D,K+1)
        return moments.astype(dtype_complex)
    def calculate_polynomial_expectation_value(self, power_string):
        S = self.S if self.S is not None else self.calculate_S()
        power = jnp.asarray(power_string, dtype=jnp.int32)           # (D,)
        moments = self.moments                                    # (n,n,D,K_MAX+1)

        idx = power[None, None, :, None]
        picked = jnp.take_along_axis(moments, idx, axis=-1)[..., 0]   # (n,n,D)
        result = jnp.prod(picked, axis=-1)                            # (n,n)
        return S * result

    def calculate_VT_poly(self, power_string):
        S = self.S if self.S is not None else self.calculate_S()
        power = jnp.asarray(power_string, dtype=jnp.int32)            # (D,)
        moments = self.moments

        k, y   = self.tilde_k, self.tilde_y
        bi, bj = self.beta_i, self.beta_j

        idx  = power[None, None, :, None]
        Mn   = jnp.take_along_axis(moments, idx, axis=-1)[..., 0]

        idx1 = jnp.clip(power - 1, 0)[None, None, :, None]
        Mn1  = jnp.take_along_axis(moments, idx1, axis=-1)[..., 0]
        Mn1  = jnp.where((power >= 1)[None, None, :], Mn1, 0)

        idx2 = jnp.clip(power - 2, 0)[None, None, :, None]
        Mn2  = jnp.take_along_axis(moments, idx2, axis=-1)[..., 0]
        Mn2  = jnp.where((power >= 2)[None, None, :], Mn2, 0)

        Mprod = jnp.prod(Mn, axis=-1)                                 # (n,n)

        inv_Mn = jnp.where(jnp.abs(Mn) > 0, 1.0/Mn, 0.0)
        r1 = Mn1 * inv_Mn
        r2 = Mn2 * inv_Mn

        nvec = power.astype(Mn.real.dtype)[None, None, :]

        base = -(2*(k**2)*(y**2) - k)
        vt_coeff = base - 2*nvec*k*y*bj*r1 - 0.5*nvec*(nvec-1)*(bj**2)*r2

        VT = S * (Mprod * jnp.sum(vt_coeff, axis=-1))
        return VT

    def calculate_TV_poly(self, power_string):
        VT = self.calculate_VT_poly(power_string)
        return VT.conj().T

    def calculate_H(self, t):
        return NotImplementedError

    def calculate_H2(self, t):
        return NotImplementedError
    
    
    def calculate_SHH2(self, t=0, params=None,splitting_type="none"):
        if params is not None:
            self.update_parameters(params)
        return self._calc_SHH2(t=t, splitting_type=splitting_type)
    #@partial(jax.jit, static_argnames=('splitting_type',), donate_argnums=(0,))
    def _calc_SHH2(self, t=0, splitting_type="none"):
        S = self.S if self.S is not None else self.calculate_S()
        if splitting_type == "none":
            H  = self.calculate_H(t)
            H2 = self.calculate_H2(t)
        else:
            H  = self.calculate_V(t)
            H2 = self.calculate_V2(t)
        return S, H, H2    
    def calculate_TV_poly(self, power_string):
        """
        Exact relation: (VT)† = TV  ⇒  TV = VT.conj().T
        """
        VT = self.calculate_VT_poly(power_string)
        return VT.conj().T
class generalPotentialSolver(ND_potentials):
    def __init__(self, params, polynomial_string):
        self.polynomial=read_string(polynomial_string)
        self.polynomial_squared=obtain_polynomial_squared(self.polynomial)
        self._poly_keys = jnp.asarray(list(self.polynomial.keys()), dtype=jnp.int32)
        self._poly_vals = jnp.asarray(list(self.polynomial.values()), dtype=dtype_complex)
        self._poly2_keys = jnp.asarray(list(self.polynomial_squared.keys()), dtype=jnp.int32)
        self._poly2_vals = jnp.asarray(list(self.polynomial_squared.values()), dtype=dtype_complex)

        K_max=4
        for key, val in self.polynomial_squared.items():
            K_max = max(K_max, max(key))
        super().__init__(params, K_max=K_max)
        self._build_cached_indices()

    def _build_cached_indices(self):
        D, Kp1 = self.D, self.moments.shape[-1]
        def _mk(keys):
            # keys: (M,D)
            dm = jnp.swapaxes(keys, 0, 1)                                  # (D,M)
            idx  = jnp.clip(dm,     0, Kp1-1)                               # (D,M)
            idx1 = jnp.clip(dm-1,   0, Kp1-1)
            idx2 = jnp.clip(dm-2,   0, Kp1-1)
            return dm, idx, idx1, idx2
        (self._dm,  self._idx,  self._idx1,  self._idx2)  = _mk(self._poly_keys)
        (self._dm2, self._i2,   self._i2_1,  self._i2_2)  = _mk(self._poly2_keys)

    def setUpIntermediates(self):
        super().setUpIntermediates()
        if hasattr(self, "_poly_keys"):
            self._build_cached_indices()
    def calculate_V(self, t=0):
        S = self.S if self.S is not None else self.calculate_S()
        M  = self._poly_vals.shape[0]
        idx = jnp.broadcast_to(self._idx,  self.moments.shape[:-1] + (M,))
        Mn  = jnp.take_along_axis(self.moments, idx, axis=-1)              # (n,n,D,M)
        Mprod = jnp.prod(Mn, axis=2)                                       # (n,n,M)
        return jnp.einsum('ijm,m->ij', S[...,None]*Mprod, self._poly_vals)

    def calculate_V2(self, t=0):
        S = self.S if self.S is not None else self.calculate_S()
        M2 = self._poly2_vals.shape[0]
        idx = jnp.broadcast_to(self._i2,  self.moments.shape[:-1] + (M2,))
        Mn  = jnp.take_along_axis(self.moments, idx, axis=-1)
        Mprod = jnp.prod(Mn, axis=2)
        return jnp.einsum('ijm,m->ij', S[...,None]*Mprod, self._poly2_vals)

    def calculate_TVVT(self, t=0):
        S = self.S if self.S is not None else self.calculate_S()
        k, y, bj = self.tilde_k, self.tilde_y, self.beta_j

        M = self._poly_vals.shape[0]
        idx  = jnp.broadcast_to(self._idx,  self.moments.shape[:-1] + (M,))
        idx1 = jnp.broadcast_to(self._idx1, self.moments.shape[:-1] + (M,))
        idx2 = jnp.broadcast_to(self._idx2, self.moments.shape[:-1] + (M,))

        Mn  = jnp.take_along_axis(self.moments, idx,  axis=-1)             # (n,n,D,M)
        Mn1 = jnp.take_along_axis(self.moments, idx1, axis=-1)
        Mn2 = jnp.take_along_axis(self.moments, idx2, axis=-1)

        # zero the invalid (n<1 / n<2) entries without branches
        Mn1 = jnp.where(self._dm[None,None,:,:] >= 1, Mn1, 0)
        Mn2 = jnp.where(self._dm[None,None,:,:] >= 2, Mn2, 0)

        inv_Mn = jnp.where(jnp.abs(Mn) > 0, 1.0/Mn, 0.0)
        r1, r2 = Mn1*inv_Mn, Mn2*inv_Mn
        nvec   = self._dm.astype(Mn.real.dtype)[None,None,:,:]

        base = -(2*(k**2)*(y**2) - k)[...,None]                            # (n,n,D,1)
        vt_coeff = base - 2*nvec*k[...,None]*y[...,None]*bj[...,None]*r1 \
                        - 0.5*nvec*(nvec-1)*(bj[...,None]**2)*r2
        term = (S[...,None] * (jnp.prod(Mn, axis=2) * jnp.sum(vt_coeff, axis=2)))  # (n,n,M)
        VT = jnp.einsum('ijm,m->ij', term, self._poly_vals)
        return VT + VT.conj().T
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
        if splitting_type=="none":
            H = self.calculate_H(t)
            H2 = self.calculate_H2(t)
        else:
            H=self.calculate_V(t)
            H2= self.calculate_V2(t)
        return S, H, H2

if __name__=="__main__":
    from jax_Rothe import *
    n=2
    D=4
    p_init = jnp.zeros((n, D, 4),dtype=dtype_real)
    p_init = p_init.at[:, :, 0].set(1/jnp.sqrt(2)) #Width parameters
    p_init = p_init.at[:, :, 0].set(1/jnp.sqrt(2)) #Width parameters
    p_init = p_init.at[0,:,2].set(2.0) #Set mu to (2,...,2)
    example_string="""
    dimension 2
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x0x0x1: 0.111803
    x1x1x1: -0.03726766666
    x0x0x0x0: 0.0007812444255625
    x1x1x1x1: 0.0007812444255625
    x0x0x1x1: 0.001562488851125
    """
    for i in range(1,n):
        p_init = p_init.at[i, :, 2].set(2+np.random.uniform(-0.5, 0.5, (D,)))          #Move to the right
    osc = generalPotentialSolver(p_init,example_string)
    S=osc.calculate_S()
    gaussian=jnp.zeros((2,D,4),dtype=dtype_real)
    gaussian=gaussian.at[0,:,0].set(1e-11) #a
    gaussian=gaussian.at[1,:,0].set(1e-11) #b
    linear_params=jnp.ones(gaussian.shape[0],dtype=dtype_complex)
    gev=calculate_Gaussian_expectation_values(p_init,gaussian,linear_params)
    print("S=",S)
    print(gev)