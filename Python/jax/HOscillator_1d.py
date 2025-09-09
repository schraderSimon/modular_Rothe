import jax.numpy as jnp
import jax
from jax.nn import softplus
from jax import jit, grad
from read_string import *
jax.config.update('jax_enable_x64', True)
import numpy as np
import jax.scipy as jsp
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.optimize import minimize
logpi = jnp.log(jnp.pi)
import numpy as np
import jax.lax as lax

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

if __name__ == "__main__":
    test_implementation_1D()
