import jax.numpy as jnp
import jax
from jax.nn import softplus
from jax import jit, grad
jax.config.update('jax_enable_x64', True)
import numpy as np
import jax.scipy as jsp
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.optimize import minimize
logpi = jnp.log(jnp.pi)
import numpy as np
class ND_potentials:
    def __init__(self, params, K_max=4):
        """
        params: (n, D, 4) -> (a, b, mu, p) per dimension
        K_max:  maximum polynomial order to precompute moments for (static under jit)
        """
        self.params = jnp.asarray(params, dtype=jnp.float64)
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
        self.params = jnp.asarray(params, dtype=jnp.float64)
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
        expo_pref = 0.5 * (D * jnp.log(jnp.pi) - jnp.sum(jnp.log(self.Gamma), axis=-1)) # (n,n)
        expo = expo_core + expo_pref  # (n,n)

        diag_expo = jnp.real(jnp.diag(expo))
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
        # Use a fixed Python int for array sizes to avoid tracers in shapes under jit
        if maximal_order is None:
            maximal_order = self.K_MAX
        K = int(maximal_order)
        n, D = self.params.shape[:2]
        moments = jnp.zeros((n, n, D, K+1), dtype=jnp.complex128)

        P, R = self.P, self.R
        M0, M1 = jnp.ones_like(P), P
        moments = moments.at[:,:,:,0].set(M0)
        moments = moments.at[:,:,:,1].set(M1)

        M_prev, M_curr = M0, M1
        for k in range(2, K+1):
            kf = jnp.asarray(k-1, dtype=P.real.dtype)
            M_next = P * M_curr + kf * R * M_prev
            M_prev, M_curr = M_curr, M_next
            moments = moments.at[:,:,:,k].set(M_curr)

        return moments

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

    def calculate_SHH2(self, t=0, params=None):
        if params is not None:
            self.update_parameters(params)
        S = self.calculate_S()
        H = self.calculate_H(t)
        H2 = self.calculate_H2(t)
        return S, H, H2
def calculate_TV_poly(self, power_string):
    """
    Exact relation: (VT)† = TV  ⇒  TV = VT.conj().T
    """
    VT = self.calculate_VT_poly(power_string)
    return VT.conj().T

class oneD_potentials():
    def __init__(self, params):
        self.params   = jnp.asarray(params,  dtype=jnp.float64)      # (n,4)
        self.n, self.N_b = self.params.shape
        assert self.N_b == 4, "each Gaussian needs 4 real parameters"
        self.mu=self.params[:,2] #Shift parameters
        self.a=self.params[:,0] #Width parameters
        self.b=self.params[:,1] #Phase parameters
        self.p=self.params[:,3] #Momentum parameters
        self.S=None
        self.setUpIntermediates()
    def setUpIntermediates(self):
        alpha = self.a**2 + 1j*self.b
        ai, aj = alpha.conj()[:, None], alpha[None, :]

        Gamma = ai + aj
        self.Gamma = Gamma
        self.exponent_contrib=-(self.p.conj()**2)[:, None]/(4*ai)-(self.p**2)[None, :]/(4*aj)
        self.tilde_k = (ai * aj) / Gamma
        mui = self.mu[:, None] - 1j * self.p[:, None] / (2 * ai)
        muj = self.mu[None, :] + 1j * self.p[None, :] / (2 * aj)
        self.tilde_y = mui - muj
        self.P= P = (ai* mui + aj* muj)/Gamma
        self.R = 1/(2 * Gamma)
        self.beta_i = ai/Gamma
        self.beta_j = aj/Gamma
        self.tysq = self.tilde_y**2
        self.tksq = self.tilde_k**2
    def update_parameters(self,params):
        self.params = jnp.asarray(params, dtype=jnp.float64)
        self.mu = self.params[:, 2]
        self.a = self.params[:, 0]
        self.b = self.params[:, 1]
        self.p = self.params[:, 3]
        self.setUpIntermediates()
    def calculate_S(self):
        S=jnp.sqrt(jnp.pi/self.Gamma)*jnp.exp(self.exponent_contrib-self.tysq*self.tilde_k)
        d = jnp.sqrt(jnp.diag(S))
        self.S = S / (d[:, None] * d[None, :])

        return self.S
    def calculate_T(self):
        if self.S is None:
            S=self.calculate_S()
        else:
            S=self.S
        T=S*(self.tilde_k-2*self.tksq*self.tysq)
        return T
    def calculate_Tsq(self):
        if self.S is None:
            S=self.calculate_S()
        else:
            S=self.S
        tk=self.tilde_k
        tksq=self.tksq
        tysq=self.tysq
        tk3=tksq*tk
        return S*(4*tksq**2*tysq**2-12*tk3*tysq+3*tksq)
    def calculate_moments(self, n):
        P, R = self.P, self.R
        M0 = jnp.ones_like(P)
        M1 = P
        if n == 0:
            return M0
        if n == 1:
            return M1

        M_prev, M_curr = M0, M1
        for k in range(1, n): # Recurrence relation 
            M_next = P * M_curr + k * R * M_prev
            M_prev, M_curr = M_curr, M_next
        return M_curr
    def calculate_x_power(self,n):
        if self.S is None:
            S=self.calculate_S()
        else:
            S=self.S
        return S*self.calculate_moments(n)
    def calculate_TVVT(self, n):
        S = getattr(self, "S", None)
        if S is None: S = self.calculate_S()
        k, y = self.tilde_k, self.tilde_y
        bi, bj = self.beta_i, self.beta_j
        Mn  = self.calculate_moments(n)
        Mn1 = self.calculate_moments(n-1) if n >= 1 else jnp.zeros_like(Mn)
        Mn2 = self.calculate_moments(n-2) if n >= 2 else jnp.zeros_like(Mn)
        term1 = -(4*k*k*y*y - 2*k) * Mn
        term2 =  2*n*k*y * (bi - bj) * Mn1
        term3 = -0.5 * n*(n-1) * (bi**2 + bj**2) * Mn2
        return S * (term1 + term2 + term3)
    def calculate_H(self,t):
        return NotImplementedError
    def calculate_H2(self,t):
        return NotImplementedError
    def calculate_SHH2(self,t=0,params=None):
        if params is not None:
            self.update_parameters(params)
        S=self.calculate_S()
        H=self.calculate_H(t)
        H2=self.calculate_H2(t)
        return S,H,H2
class HOscillator(oneD_potentials):
    def __init__(self, params, omega=1.0, name="HOsc_1D"):
        super().__init__(params)
        self.omega = omega
    def calculate_H(self,t=0):
        return self.calculate_T()+0.5*self.calculate_x_power(2)*self.omega**2 #Harmonic oscillator Hamiltonian
    def calculate_H2(self,t=0):
        return self.calculate_Tsq()+0.25*self.calculate_x_power(4)*self.omega**4 +0.5*self.calculate_TVVT(2)*self.omega**2
class HOscillator_ND(ND_potentials):
    def __init__(self, params, omega=1.0, name="HOsc_ND", K_max=4):
        super().__init__(params, K_max=K_max)
        self.omega = omega

    def calculate_H(self, t=0):
        D = self.D
        H = self.calculate_T()
        for i in range(D):
            polynom = jnp.asarray([0]*i + [2] + [0]*(D-i-1))
            H += 0.5 * self.calculate_polynomial_expectation_value(polynom) * self.omega**2
        return H

    def calculate_H2(self, t=0):
        D = self.D
        H2 = self.calculate_Tsq()
        for i in range(D):
            for j in range(D):
                polynom = jnp.zeros(D, dtype=int)
                polynom = polynom.at[i].add(2).at[j].add(2)
                polynom = jnp.asarray(polynom)
                H2 += 0.25 * self.calculate_polynomial_expectation_value(polynom) * self.omega**4

        TVVT = jnp.zeros_like(H2)
        for i in range(D):
            polynom = jnp.asarray([0]*i + [2] + [0]*(D-i-1))
            TVVT += 0.5 * self.calculate_VT_poly(polynom) * self.omega**2
        TVVT = TVVT + jnp.conj(TVVT.T)

        return H2 + TVVT
def eigh_canonical(S, H, eps=1e-12):
    S = 0.5*(S + S.conj().T)
    H = 0.5*(H + H.conj().T)
    w, U = jnp.linalg.eigh(S)
    wmax = jnp.max(jnp.real(w))
    keep = jnp.real(w) > eps*wmax
    U1, w1 = U[:, keep], w[keep]
    X = U1 * (1.0 / jnp.sqrt(w1))[None, :]
    H_orth = X.conj().T @ H @ X
    e, V = jnp.linalg.eigh(H_orth)
    C = X @ V
    return e, C
def _symh(A):  # Hermitianize (numerical hygiene)
    return 0.5 * (A + A.conj().T)

def groundstate_energy(params, omega=1.0, eps=1e-12):
    ho = HOscillator(params, omega=omega)
    S,H,H2=ho.calculate_SHH2()
    S = _symh(S)
    H = _symh(H)
    # 1) diagonal pre-normalization: D^{-1/2} S D^{-1/2} has unit diagonal
    d = jnp.clip(jnp.real(jnp.diag(S)), a_min=eps, a_max=None)
    Dm12 = 1.0 / jnp.sqrt(d)
    S1 = (Dm12[:, None] * S) * Dm12[None, :]
    H1 = (Dm12[:, None] * H) * Dm12[None, :]
    # 2) exact canonical orthonormalization: S1^{-1/2}
    w, U = jnp.linalg.eigh(S1)
    w = jnp.clip(jnp.real(w), a_min=eps, a_max=None)
    X = U * (1.0 / jnp.sqrt(w))[None, :]           # S1^{-1/2}
    H_orth = X.conj().T @ H1 @ X                   # Hermitian
    e = jnp.linalg.eigvalsh(H_orth)                # sorted asc
    return jnp.real(e[0])
# JIT + grad on the packed parameter vector
def _pack(P): return jnp.ravel(P)
def _unpack(theta, n): return theta.reshape((n,4))

def make_energy_funcs(n, omega=1.0, eps=1e-12):
    f_theta = lambda theta: groundstate_energy(_unpack(theta, n), omega=omega, eps=eps)
    f_jit = jit(f_theta)
    g_jit = jit(grad(f_theta))
    def f_np(theta): return np.asarray(f_jit(jnp.asarray(theta))).item()
    def g_np(theta): return np.asarray(g_jit(jnp.asarray(theta)))
    return f_np, g_np

def find_groundstate(params0, omega=1.0, eps=1e-12):
    params0 = jnp.asarray(params0, dtype=jnp.float64)
    n = params0.shape[0]
    f_np, g_np = make_energy_funcs(n, omega=omega, eps=eps)
    x0 = np.asarray(_pack(params0))
    res = minimize(f_np, x0, jac=g_np, method="BFGS",
                   options=dict(maxiter=500))
    return res.x.reshape(n,4), res.fun, res
def groundstate_variance(params, omega=1.0, eps=1e-12):
    ho   = HOscillator(params, omega=omega)
    S,H,H2=ho.calculate_SHH2()
    # Canonical orthonormalization: S1^{-1/2}
    w, U = jnp.linalg.eigh(S)
    w    = jnp.clip(jnp.real(w), a_min=eps, a_max=None)
    X    = U * (1.0 / jnp.sqrt(w))[None, :]

    H_orth   = X.conj().T @ H @ X
    H2_orth = X.conj().T @ H2 @ X

    # Lowest eigenpair of H_orth
    e, V = jnp.linalg.eigh(H_orth)
    v0   = V[:, 0]

    # Var = <H^2> - <H>^2, in the orthonormalized space
    var = jnp.real(v0.conj() @ (H2_orth @ v0)) - (jnp.real(e[0])**2)
    return var
def make_variance_funcs(n, omega=1.0, eps=1e-12):
    f_theta = lambda theta: groundstate_variance(_unpack(theta, n), omega=omega, eps=eps)
    f_jit   = jit(f_theta)
    g_jit   = jit(grad(f_theta))
    def f_np(theta): return np.asarray(f_jit(jnp.asarray(theta))).item()
    def g_np(theta): return np.asarray(g_jit(jnp.asarray(theta)))
    return f_np, g_np

def minimize_variance(params0, omega=1.0, eps=1e-12):
    params0 = jnp.asarray(params0, dtype=jnp.float64)
    n = params0.shape[0]
    f_np, g_np = make_variance_funcs(n, omega=omega, eps=eps)
    x0 = np.asarray(_pack(params0))
    # bounds: a_i > 0, others free
    res = minimize(f_np, x0, jac=g_np, method="BFGS",
                   options=dict(maxiter=500))
    return res.x.reshape(n, 4), res.fun, res
def test_implementation_1D():
    params=jnp.array([[1/jnp.sqrt(2), 0.0, 0.0, 0.0]]) # Ground state of HO
    hoscillator=HOscillator(params)
    S,H,H2=hoscillator.calculate_SHH2() 
    energy=H[0,0]
    assert jnp.isclose(energy, 0.5)
    variance=H2[0,0]-energy**2
    assert jnp.isclose(variance, 0.0)
    num=10
    params=jnp.zeros((num,4))
    params = params.at[:, 0].set(1/jnp.sqrt(3)) #Width parameters
    params = params.at[:num//2,  3].set(jnp.linspace(-1,1,num//2)) #First dimension gets different mu parameter
    params = params.at[:num//2,  4].set(1) #Shift in momentum
    params = params.at[:num//2,  1].set(0.1)
    params = params.at[num//2:,  3].set(jnp.linspace(1,-1,num//2))
    params = params.at[num//2:,  4].set(-1)
    params = params.at[num//2:,  1].set(-0.1)
    hoscillator=HOscillator(params)
    S,H,H2=hoscillator.calculate_SHH2()
    eigvals,eigvecs=eigh_canonical(S,H)
    print(eigvals)
    """
    H2_diagonalbasis=eigvecs.conj().T @ H2 @ eigvecs
    variances=jnp.diag(H2_diagonalbasis-jnp.diag(eigvals)**2)
    assert jnp.isclose(variances[0],0)

    assert jnp.isclose(eigvals[0],0.5)
    startguess=jnp.array([[1/jnp.sqrt(3), 2.0, 1.0, 1.0],[1, 0.0, -1.0, 1.0]])


    params,energy,res=find_groundstate(startguess)
    emin=groundstate_energy(params)
    assert jnp.isclose(emin, 0.5)
    print("Energy of e-minimized state:", emin)
    print("Variance of e-minimized state:", groundstate_variance(params))

    P_opt, var_opt, res = minimize_variance(startguess, omega=1.0)
    print("variance* =", var_opt)
    print("Energy of variance-minimized state:", groundstate_energy(P_opt))
    assert jnp.isclose(var_opt, 0)
    """
def test_implementation_ND(D=3):
    params=jnp.zeros((1, D, 4))
    params = params.at[0, :, 0].set(1/jnp.sqrt(2)) #Width parameters
    hoscillator=HOscillator_ND(params)
    S,H,H2=hoscillator.calculate_SHH2()
    energy=H[0,0]
    kinetic_energy=hoscillator.calculate_T()[0,0]
    variance=H2[0,0]-energy**2
    print(kinetic_energy,energy,H2[0,0])
    assert jnp.isclose(energy, 0.5*D)
    assert jnp.isclose(variance, 0.0)
    num=10
    params=jnp.zeros((num,D,4))
    params = params.at[:, :, 0].set(1/jnp.sqrt(2)) #Width parameters
    params = params.at[:, 0, 0].set(1/jnp.sqrt(3))
    params = params.at[:num//2, 0, 3].set(jnp.linspace(-1,1,num//2)) #First dimension gets different mu parameter
    params = params.at[:num//2, 0, 4].set(1) #Shift in momentum
    params = params.at[:num//2, 0, 1].set(0.1)
    params = params.at[num//2:, 0, 3].set(jnp.linspace(1,-1,num//2))
    params = params.at[num//2:, 0, 4].set(-1)
    params = params.at[num//2:, 0, 1].set(-0.1)
    hoscillator=HOscillator_ND(params)
    S,H,H2=hoscillator.calculate_SHH2()
    eigvals,eigvecs=eigh_canonical(S,H)
    H2_diagonalbasis=eigvecs.conj().T @ H2 @ eigvecs
    variances=jnp.diag(H2_diagonalbasis-jnp.diag(eigvals)**2)
    print(eigvals)
    print(variances)
    assert jnp.isclose(variances[0],0)

    assert jnp.isclose(eigvals[0],0.5*D)
    
if __name__ == "__main__":
    test_implementation_1D()
    test_implementation_ND(10)