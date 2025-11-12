from abc import ABC, abstractmethod
from typing import cast

import jax
import jax.numpy as jnp
from jax import grad, jit

jax.config.update("jax_enable_x64", True)

import numpy as np
from scipy.optimize import minimize

logpi = jnp.log(jnp.pi)
import numpy as np

from utils import _symh, divide_outer


class oneD_potentials(ABC):
    def __init__(self, params):
        self.params = jnp.asarray(params, dtype=jnp.float64)  # (n,4)
        self.n, self.N_b = self.params.shape
        assert self.N_b == 4, "each Gaussian needs 4 real parameters"
        self.mu = self.params[:, 2]  # Shift parameters
        self.a = self.params[:, 0]  # Width parameters
        self.b = self.params[:, 1]  # Phase parameters
        self.p = self.params[:, 3]  # Momentum parameters
        self.S = None
        self.setUpIntermediates()

    def setUpIntermediates(self):
        alpha = self.a**2 + 1j * self.b

        alpha_conj = alpha.conj()
        Gamma = jnp.add.outer(alpha_conj, alpha)
        self.tilde_k = jnp.multiply.outer(alpha_conj, alpha) / Gamma

        muj = self.mu + 1j * self.p / (2 * alpha)
        mui = muj.conj()
        self.tilde_y = jnp.subtract.outer(mui, muj)
        self.P = jnp.add.outer(alpha_conj * mui, alpha * muj) / Gamma
        self.R = 1 / (2 * Gamma)
        self.tysq = self.tilde_y**2
        self.tksq = self.tilde_k**2
        ones = jnp.ones_like(alpha)

        self.beta_j = jnp.outer(alpha, ones) / Gamma
        self.beta_i = self.beta_j.conj()
        self.Gamma = Gamma

        exponent = self.p**2 / (4 * alpha)
        self.exponent_contrib = -jnp.add.outer(exponent.conj(), exponent)

    def update_parameters(self, params):
        self.params = jnp.asarray(params, dtype=jnp.float64)
        self.mu = self.params[:, 2]
        self.a = self.params[:, 0]
        self.b = self.params[:, 1]
        self.p = self.params[:, 3]
        self.setUpIntermediates()

    def calculate_S(self):
        S = jnp.sqrt(jnp.pi / self.Gamma) * jnp.exp(
            self.exponent_contrib - self.tysq * self.tilde_k
        )
        d = jnp.sqrt(jnp.diag(S))  # Normalization constants

        # Normalize
        self.S = S / jnp.multiply.outer(d, d)
        return self.S

    def calculate_T(self):
        if self.S is None:
            S = self.calculate_S()
        else:
            S = self.S
        T = S * (self.tilde_k - 2 * self.tksq * self.tysq)
        return T

    def calculate_Tsq(self):
        if self.S is None:
            S = self.calculate_S()
        else:
            S = self.S
        tk = self.tilde_k
        tksq = self.tksq
        tysq = self.tysq
        tk3 = tksq * tk
        return S * (4 * tksq**2 * tysq**2 - 12 * tk3 * tysq + 3 * tksq)

    def calculate_moments(self, n):
        P, R = self.P, self.R
        M0 = jnp.ones_like(P)
        M1 = P
        if n == 0:
            return M0
        if n == 1:
            return M1

        M_prev, M_curr = M0, M1
        for k in range(1, n):  # Recurrence relation
            M_next = P * M_curr + k * R * M_prev
            M_prev, M_curr = M_curr, M_next
        return M_curr

    def calculate_x_power(self, n):
        if self.S is None:
            S = self.calculate_S()
        else:
            S = self.S
        return S * self.calculate_moments(n)

    def calculate_TVVT(self, n):
        S = getattr(self, "S", None)
        if S is None:
            S = self.calculate_S()
        k, y = self.tilde_k, self.tilde_y
        bi, bj = self.beta_i, self.beta_j
        Mn = self.calculate_moments(n)
        Mn1 = self.calculate_moments(n - 1) if n >= 1 else jnp.zeros_like(Mn)
        Mn2 = self.calculate_moments(n - 2) if n >= 2 else jnp.zeros_like(Mn)
        term1 = -(4 * k * k * y * y - 2 * k) * Mn
        term2 = 2 * n * k * y * (bi - bj) * Mn1
        term3 = -0.5 * n * (n - 1) * (bi**2 + bj**2) * Mn2
        return S * (term1 + term2 + term3)

    @abstractmethod
    def calculate_H(self, t):
        raise NotImplementedError

    @abstractmethod
    def calculate_H2(self, t):
        raise NotImplementedError

    def calculate_SHH2(self, t=0, params=None):
        if params is not None:
            self.update_parameters(params)
        S = self.calculate_S()
        H = self.calculate_H(t)
        H2 = self.calculate_H2(t)
        return S, H, H2


class HOscillator(oneD_potentials):
    def __init__(self, params, omega=1.0, name="HOsc_1D"):
        super().__init__(params)
        self.omega = omega

    def calculate_H(self, t=0):
        return (
            self.calculate_T() + 0.5 * self.calculate_x_power(2) * self.omega**2
        )  # Harmonic oscillator Hamiltonian

    def calculate_H2(self, t=0):
        return (
            self.calculate_Tsq()
            + 0.25 * self.calculate_x_power(4) * self.omega**4
            + 0.5 * self.calculate_TVVT(2) * self.omega**2
        )


def _unpack(theta: jnp.ndarray, n: int) -> jnp.ndarray:
    return theta.reshape((n, 4))


def make_energy_funcs(n, omega=1.0, eps=1e-12):
    f_theta = lambda theta: groundstate_energy(_unpack(theta, n), omega=omega, eps=eps)
    f_jit = jit(f_theta)
    g_jit = jit(grad(f_theta))

    def f_np(theta):
        return np.asarray(f_jit(jnp.asarray(theta))).item()

    def g_np(theta):
        return np.asarray(g_jit(jnp.asarray(theta)))

    return f_np, g_np


def find_groundstate(params0, omega=1.0, eps=1e-12):
    params0 = jnp.asarray(params0, dtype=jnp.float64)
    n = params0.shape[0]
    f_np, g_np = make_energy_funcs(n, omega=omega, eps=eps)
    x0 = np.asarray(jnp.ravel(params0))
    res = minimize(f_np, x0, jac=g_np, method="BFGS", options=dict(maxiter=500))
    return res.x.reshape(n, 4), res.fun, res


def canonical_orthogonalization(S, H, eps=1e-12):
    pass
    w, U = jnp.linalg.eigh(S)
    w = jnp.clip(jnp.real(w), a_min=eps, a_max=None)
    X = U * (1.0 / jnp.sqrt(w))[None, :]
    return X


def calculate_variance(state: jnp.ndarray, H2: jnp.ndarray, e: float) -> jnp.ndarray:
    var = jnp.real(state.conj() @ (H2 @ state)) - (jnp.real(e) ** 2)
    return var


def groundstate_variance(params, omega=1.0, eps=1e-12):
    ho = HOscillator(params, omega=omega)
    S, H, H2 = ho.calculate_SHH2()

    H_orth = canonical_orthogonalization(S, H, eps)
    H2_orth = canonical_orthogonalization(S, H2, eps)

    # Lowest eigenpair of H_orth
    e, V = jnp.linalg.eigh(H_orth)
    v0 = V[:, 0]

    # Var = <H^2> - <H>^2, in the orthonormalized space
    var = calculate_variance(v0, H2_orth, e[0])
    return var


def make_variance_funcs(n, omega=1.0, eps=1e-12):
    f_theta = lambda theta: groundstate_variance(_unpack(theta, n), omega=omega, eps=eps)
    f_jit = jit(f_theta)
    g_jit = jit(grad(f_theta))

    def f_np(theta):
        return np.asarray(f_jit(jnp.asarray(theta))).item()

    def g_np(theta):
        return np.asarray(g_jit(jnp.asarray(theta)))

    return f_np, g_np


def minimize_variance(params0, omega=1.0, eps=1e-12):
    params0 = jnp.asarray(params0, dtype=jnp.float64)
    n = params0.shape[0]
    f_np, g_np = make_variance_funcs(n, omega=omega, eps=eps)
    x0 = np.asarray(jnp.ravel(params0))
    # bounds: a_i > 0, others free
    res = minimize(f_np, x0, jac=g_np, method="BFGS", options=dict(maxiter=500))
    return res.x.reshape(n, 4), res.fun, res


def groundstate_energy(params, omega=1.0, eps=1e-12):
    ho = HOscillator(params, omega=omega)
    S, H, H2 = ho.calculate_SHH2()
    S = _symh(S)
    H = _symh(H)
    # 1) diagonal pre-normalization: D^{-1/2} S D^{-1/2} has unit diagonal
    d = jnp.clip(jnp.real(jnp.diag(S)), a_min=eps, a_max=None)
    Dm12 = 1.0 / jnp.sqrt(d)
    Dm_outer = jnp.outer(Dm12, Dm12)
    S1 = S * Dm_outer
    H1 = H * Dm_outer

    # 2) exact canonical orthonormalization: S1^{-1/2}
    w, U = jnp.linalg.eigh(S1)
    w = jnp.clip(jnp.real(w), a_min=eps, a_max=None)
    X = U * jnp.outer(jnp.ones(U.shape[0]), 1.0 / jnp.sqrt(w))  # S1^{-1/2}
    H_orth = X.conj().T @ H1 @ X  # Hermitian
    e = jnp.linalg.eigvalsh(H_orth)  # sorted asc
    return jnp.real(e[0])
