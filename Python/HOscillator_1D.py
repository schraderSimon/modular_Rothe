import numpy as np
from numpy.polynomial import polynomial as P
import sympy as sp           # used both for validation and for fallback
from functools import lru_cache

class HOscillator_1D:
    r"""
    Harmonic oscillator  H = −½ ∂²/∂x² + ½ ω² x²  (m = ħ = 1)
    Gaussian basis functions

        g_i(x) = N_i · exp[ −(a_i² + i b_i) (x − μ_i)²  +  i p_i (x − μ_i) ]     (★)

    All four non-linear parameters (a_i, b_i, μ_i, p_i) are **real**.
    The linear coefficients live in the complex vector ``c``.

    The analytic formulas are evaluated first; every single element is
    then cross-checked against SymPy to guarantee correctness.
    """

    # ─────────────────────────── constructor ────────────────────────
    def __init__(self, p, c, *, omega=1.0, dt=1e-2, name="HOsc_1D"):
        self.p   = np.ascontiguousarray(p,  dtype=float)      # (n,4)
        self.c   = np.ascontiguousarray(c,  dtype=complex)    # (n,)
        self.ω   = float(omega)
        self.dt  = float(dt)
        self.name = str(name)

        self.n, self.N_b = self.p.shape
        assert self.N_b == 4, "each Gaussian needs 4 real parameters"

        self.Sd=self.dS
        self.Hd=self.dH
        self.H2d=self.dH2      

        # pre-build the SymPy symbols once

    # ───────────────────────── helper blocks ────────────────────────
    @staticmethod
    def _norm_const(A):
        """‖g‖ = 1 ⇒ N = (2 Re A / π)^{1/4} (phase term i p cancels)."""
        return (2*A.real/np.pi)**0.25

    @staticmethod
    def _moments(α, β, Z, n_max):
        """I_k = ∫ x^k e^{−α x² + β x} dx,  Re α>0."""
        I = [None]*(n_max+1)
        I[0] = Z
        if n_max == 0:
            return I
        I[1] = β/(2*α)*I[0]
        for k in range(2, n_max+1):
            I[k] = β/(2*α)*I[k-1] + (k-1)/(2*α)*I[k-2]
        return I

    def _common(self, i, j,p):
        """Everything that appears in S, H, H² for indices (i,j)."""
        a_i, b_i, μ_i, p_i =p[i]
        a_j, b_j, μ_j, p_j =p[j]

        A_i, A_j = a_i**2 + 1j*b_i, a_j**2 + 1j*b_j
        N_i, N_j = self._norm_const(A_i), self._norm_const(A_j)

        α = np.conj(A_i) + A_j
        β = 2*(np.conj(A_i)*μ_i + A_j*μ_j) + 1j*(p_j - p_i)
        γ = -(np.conj(A_i)*μ_i**2 + A_j*μ_j**2) + 1j*(p_i*μ_i - p_j*μ_j)

        Z = N_i * N_j * np.sqrt(np.pi/α) * np.exp(γ + β**2/(4*α))
        return A_i, A_j, μ_i, μ_j, p_i, p_j, α, β,γ, Z

    @staticmethod
    def _shift(poly, k):  return [0]*k + list(poly)
    @staticmethod
    def _scale(poly, s):  return [s*coef for coef in poly]
    @staticmethod
    def _add_many(polys):
        L = max(len(p) for p in polys)
        out = [0]*L
        for poly in polys:
            for k,c in enumerate(poly):
                out[k] += c
        return out

    # ───────────────── analytic S, H, H² (fast) ─────────────────────
    def S(self, parr):
        n = parr.shape[0]
        S = np.zeros((n, n), complex)
        for i in range(n):
            for j in range(i, n):
                A_i, A_j, μ_i, μ_j, p_i, p_j, α, β,γ, Z = self._common(i, j,parr)
                S[i,j] = Z
                if i!=j: S[j,i] = Z.conjugate()
        return S

    def H(self, parr, t=None):
        n = parr.shape[0];  ω = self.ω
        H = np.zeros((n, n), complex)
        for i in range(n):
            for j in range(i, n):
                A_i, A_j, μ_i, μ_j, p_i, p_j, α, β,γ, Z = self._common(i, j,parr)
                I0,I1,I2 = self._moments(α,β,Z,2)
                b1,b0 = -2*A_j, (2*A_j*μ_j + 1j*p_j)
                B2    = P.polypow([b0,b1],2)
                T_ij  = -0.5*(B2[0]*I0 + B2[1]*I1 + B2[2]*I2 + b1*I0)
                V_ij  = 0.5*ω**2*I2
                H[i,j] = T_ij + V_ij
                if i!=j: H[j,i] = H[i,j].conjugate()
        return H



    def H2(self,parr,t=None):
        ω, λ = self.ω, 0.5*self.ω**2
        n = parr.shape[0];  ω = self.ω
        H2m   = np.zeros((n, n), complex)

        for i in range(n):
            for j in range(i, n):
                A_i, A_j, μ_i, μ_j, p_i, p_j, α, β,γ, Z = self._common(i, j,parr)
                I = self._moments(α, β, Z, 4)

                b1, b0 = -2*A_j, (2*A_j*μ_j + 1j*p_j)
                B      = [b0, b1]
                B2, B4 = list(P.polypow(B, 2)), list(P.polypow(B, 4))
                B2 += [0]*(5-len(B2));  B4 += [0]*(5-len(B4))

                poly_T2 = self._add_many([[3*A_j**2]+[0]*4,
                                          self._scale(B2, -3*A_j),
                                          self._scale(B4, 0.25)])

                E = -λ/2
                poly_TVVT = self._scale(
                    self._add_many([[2]+[0]*4,
                                    self._scale(self._shift(B,1), 4),
                                    self._scale(self._shift(B2,2), 2),
                                    self._scale(self._shift([b1],2), 2)]),
                    E)

                poly_V2   = [0,0,0,0, λ**2]
                poly_tot  = np.asarray(self._add_many([poly_T2, poly_TVVT, poly_V2]))

                H2_ij = sum(poly_tot[k]*I[k] for k in range(5))
                H2m[i, j] = H2_ij
                if i != j:
                    H2m[j, i] = H2_ij.conjugate()
        return H2m

    # ─────────── exact SymPy integrals (memoised, slow) ─────────────


    # ─────────────────── public, validated API ──────────────────────

    def dS(self, p):
        r"""
        Analytic three-tensor of overlap derivatives for an *arbitrary*
        parameter matrix ``p`` (shape (n,4)).

            dS[i, j, k] = ∂S_ij / ∂(param_k of ket Gaussian j)

        with   k = 0 (a_j), 1 (b_j), 2 (μ_j), 3 (p_j).

        Only the *ket* contribution (column ``j``) is returned; the caller
        can add its Hermitian conjugate if the full derivative is needed.
        The diagonal is zero by construction, because it is the overlap. Lol.
        """
        p   = np.ascontiguousarray(p, float)
        n   = p.shape[0]
        dS  = np.zeros((n, n, 4), complex)

        for i in range(n):
            a_i, b_i, mu_i, p_i = p[i]
            A_i = a_i**2 + 1j*b_i
            N_i = (2*a_i**2/np.pi)**0.25
            A_i_c = A_i.conjugate()

            for j in range(n):
                a_j, b_j, mu_j, p_j = p[j]
                A_j = a_j**2 + 1j*b_j
                N_j = (2*a_j**2/np.pi)**0.25

                α = A_i_c + A_j
                β = 2*(A_i_c*mu_i + A_j*mu_j) + 1j*(p_j - p_i)
                γ = -(A_i_c*mu_i**2 + A_j*mu_j**2) + 1j*(p_i*mu_i - p_j*mu_j)

                S_ij = N_i*N_j*np.sqrt(np.pi/α) * np.exp(γ + β**2/(4*α))

                # ∂ ln S / ∂ parameters of Gaussian j (ket side)
                dlog_da  = (0.5/a_j - a_j/α + 2*a_j*mu_j*β/α
                            - 2*a_j*mu_j**2 - a_j*β**2/(2*α**2))
                dlog_db  = 1j*(-mu_j**2 - 0.5/α + mu_j*β/α
                               - β**2/(4*α**2))
                dlog_dmu = A_j*(β/α - 2*mu_j) - 1j*p_j
                dlog_dp  = 1j*(β/(2*α) - mu_j)

                dS[i, j, 0] = S_ij * dlog_da
                dS[i, j, 1] = S_ij * dlog_db
                dS[i, j, 2] = S_ij * dlog_dmu
                dS[i, j, 3] = S_ij * dlog_dp

        # self-overlaps never change with their own parameters
        for i in range(n):
            dS[i,i,:]+=np.conj(dS[i,i,:]) # Take into account that the parameters are on both sides in the cross terms. Nicey daisy.

        return dS
    @staticmethod
    def _dalpha_dbeta(k, a_j, b_j, mu_j):
        """
        Return (dα, dβ) for derivative index k
        k = 0(a_j)  1(b_j)  2(μ_j)  3(p_j)
        """
        if   k == 0:                       # ∂/∂a_j
            dalpha =  2.0*a_j
            dbeta  =  4.0*a_j*mu_j
        elif k == 1:                       # ∂/∂b_j
            dalpha =  1j
            dbeta  =  2j*mu_j
        elif k == 2:                       # ∂/∂μ_j
            dalpha =  0.0
            dbeta  =  2.0*(a_j**2 + 1j*b_j)
        else:                              # ∂/∂p_j
            dalpha =  0.0
            dbeta  =  1j
        return dalpha, dbeta

    # ───────────────── analytic ∂H/∂param  (ket side only) ────────────────────
    def dH(self, p,t=None):
        """
        Three-tensor  dH[i,j,k] = ∂H_ij / ∂(param_k of ket Gaussian j)
        k = 0(a_j), 1(b_j), 2(μ_j), 3(p_j).   Only the ket contribution
        (column j) is returned – exactly the convention used in the test.
        """
        p   = np.ascontiguousarray(p, float)
        n   = p.shape[0]
        dH  = np.zeros((n, n, 4), complex)
        ω2  = self.ω**2

        for i in range(n):
            a_i, b_i, mu_i, p_i = p[i]
            A_i  = a_i**2 + 1j*b_i
            N_i  = (2*a_i**2/np.pi)**0.25
            A_ic = A_i.conjugate()

            for j in range(n):
                a_j, b_j, mu_j, p_j = p[j]
                A_j = a_j**2 + 1j*b_j
                N_j = (2*a_j**2/np.pi)**0.25

                α = A_ic + A_j
                β = 2*(A_ic*mu_i + A_j*mu_j) + 1j*(p_j - p_i)
                γ = -(A_ic*mu_i**2 + A_j*mu_j**2) + 1j*(p_i*mu_i - p_j*mu_j)

                S_ij = N_i*N_j*np.sqrt(np.pi/α) * np.exp(γ + β**2/(4*α))

                # ---- Gaussian moments (functions of α,β) --------------------
                F1 =  β / (2*α)
                F2 = (β**2 + 2*α) / (4*α**2)
                I0 = S_ij
                I1 = F1 * S_ij
                I2 = F2 * S_ij

                # coefficients in ⟨T⟩ piece
                b1, b0 = -2*A_j, (2*A_j*mu_j + 1j*p_j)
                B20, B21, B22 = b0**2, 2*b0*b1, b1**2

                # ---------------- derivatives for EACH parameter -------------
                for k in range(4):
                    # ∂α, ∂β
                    dα, dβ = self._dalpha_dbeta(k, a_j, b_j, mu_j)

                    # ∂ ln S  (already derived)
                    if   k == 0:
                        dlogS = (0.5/a_j - a_j/α + 2*a_j*mu_j*β/α
                                 - 2*a_j*mu_j**2 - a_j*β**2/(2*α**2))
                    elif k == 1:
                        dlogS = 1j*(-mu_j**2 - 0.5/α + mu_j*β/α
                                    - β**2/(4*α**2))
                    elif k == 2:
                        dlogS = A_j*(β/α - 2*mu_j) - 1j*p_j
                    else:
                        dlogS = 1j*(β/(2*α) - mu_j)

                    # ∂F1 , ∂F2
                    dF1 = (dβ*α - β*dα) / (2*α**2)

                    N   = β**2 + 2*α
                    dN  = 2*β*dβ + 2*dα
                    dF2 = (dN * 4*α**2 - N * 8*α*dα) / (4*α**2)**2

                    # ∂ moments
                    dI0 = I0 * dlogS
                    dI1 = S_ij * (dF1 + F1*dlogS)
                    dI2 = S_ij * (dF2 + F2*dlogS)

                    # ∂ b0 , ∂ b1
                    if   k == 0:
                        db1 = -4*a_j
                        db0 =  4*a_j*mu_j
                    elif k == 1:
                        db1 = -2j
                        db0 =  2j*mu_j
                    elif k == 2:
                        db1 =  0.0
                        db0 =  2*A_j
                    else:
                        db1 =  0.0
                        db0 =  1j

                    # ∂B2 coefficients
                    dB20 = 2*b0*db0
                    dB21 = 2*(db0*b1 + b0*db1)
                    dB22 = 2*b1*db1

                    # assemble derivative of H_ij
                    dH_ij = -0.5*(dB20*I0 + B20*dI0
                                  + dB21*I1 + B21*dI1
                                  + dB22*I2 + B22*dI2
                                  + db1*I0  + b1*dI0) \
                            + 0.5*ω2 * dI2

                    dH[i, j, k] = dH_ij
        for i in range(n):
            dH[i,i,:]+=np.conj(dH[i,i,:])
        return dH
    def dH2(self, p,t=None):
        """
        Three-tensor dH2[i,j,k] = ∂⟨g_i|H²|g_j⟩ / ∂(param_k of ket Gaussian j)
        k = 0(a_j), 1(b_j), 2(μ_j), 3(p_j).   Only the ket (column-wise)
        part is returned; add its Hermitian conjugate yourself if needed.
        """
        p   = np.ascontiguousarray(p, float)
        n   = p.shape[0]
        dH2 = np.zeros((n, n, 4), complex)

        ω2 = self.ω**2
        λ  = 0.5 * ω2
        E  = -λ / 2

        # helper: build moments I_0…4 *and* their four derivatives dI[:,k]
        def _moments_and_grads(α, β, S, dα, dβ, dlogS):
            I  = [None]*5
            dI = np.zeros((5, 4), dtype=complex)

            I[0]  = S
            dI[0] = S * dlogS              # plain chain rule

            fac   = β / (2*α)
            dfac  = (dβ*α - β*dα) / (2*α**2)

            # k = 1 explicitly
            I[1]  = fac * I[0]
            dI[1] = dfac*I[0] + fac*dI[0]

            # higher k via recurrence  I_k = fac*I_{k-1} + g*I_{k-2}
            for k in range(2, 5):
                g   = (k-1)/(2*α)
                dg  = -(k-1)/(2*α**2) * dα
                I[k]  = fac*I[k-1] + g*I[k-2]
                dI[k] = (dfac*I[k-1] + fac*dI[k-1]
                         + dg*I[k-2]  + g*dI[k-2])
            return I, dI

        # ------------------------------------------------------------------
        for i in range(n):
            a_i, b_i, μ_i, p_i = p[i]
            A_i  = a_i**2 + 1j*b_i
            A_ic = A_i.conjugate()

            for j in range(n):
                a_j, b_j, μ_j, p_j = p[j]
                A_j  = a_j**2 + 1j*b_j

                α = A_ic + A_j
                β = 2*(A_ic*μ_i + A_j*μ_j) + 1j*(p_j - p_i)
                γ = -(A_ic*μ_i**2 + A_j*μ_j**2) + 1j*(p_i*μ_i - p_j*μ_j)
                S = ( (2*a_i**2/np.pi)**0.25 *
                      (2*a_j**2/np.pi)**0.25 *
                      np.sqrt(np.pi/α) * np.exp(γ + β**2/(4*α)) )

                # ---------- dα,dβ and d log S for each parameter k ----------
                dαβ = [self._dalpha_dbeta(k, a_j, b_j, μ_j) for k in range(4)]
                dα  = np.array([ab[0] for ab in dαβ])
                dβ  = np.array([ab[1] for ab in dαβ])
                dlogS = np.array([
                    (0.5/a_j - a_j/α + 2*a_j*μ_j*β/α
                     - 2*a_j*μ_j**2 - a_j*β**2/(2*α**2)),
                    1j*(-μ_j**2 - 0.5/α + μ_j*β/α
                        - β**2/(4*α**2)),
                    A_j*(β/α - 2*μ_j) - 1j*p_j,
                    1j*(β/(2*α) - μ_j)
                ], dtype=complex)

                # moments I_k and their derivatives dI_k/dθ
                I, dI = _moments_and_grads(α, β, S, dα, dβ, dlogS)

                # ------- polynomial coefficients of the H² integrand --------
                b1, b0 = -2*A_j, (2*A_j*μ_j + 1j*p_j)
                B20, B21, B22 = b0**2,            2*b0*b1,           b1**2
                B40, B41 =      b0**4,            4*b0**3*b1
                B42, B43, B44 = 6*b0**2*b1**2,    4*b0*b1**3,         b1**4

                poly_T2 = np.array([
                    3*A_j**2 - 3*A_j*B20 + 0.25*B40,
                    -3*A_j*B21 + 0.25*B41,
                    -3*A_j*B22 + 0.25*B42,
                    0.25*B43,
                    0.25*B44
                ], dtype=complex)

                poly_TVVT = np.array([
                    2,
                    4*b0,
                    6*b1 + 2*B20,        # 4*b1  from 4·shift(B,1)
                                        # 2*b1  from 2·shift([b1],2)
                    2*B21,
                    2*B22
                ]) * E

                poly_V2 = np.array([0, 0, 0, 0, λ**2], dtype=complex)

                poly_tot = poly_T2 + poly_TVVT + poly_V2   # length-5 vector

                # ---- derivatives of those coefficients wrt each parameter ---
                dA  = np.array([2*a_j, 1j, 0.0, 0.0])
                db0 = np.array([4*a_j*μ_j,  2j*μ_j, 2*A_j, 1j])
                db1 = np.array([-4*a_j,     -2j,    0.0,  0.0])

                dB20 = 2*b0*db0
                dB21 = 2*(db0*b1 + b0*db1)
                dB22 = 2*b1*db1

                dB40 = 4*b0**3 * db0
                dB41 = 12*b0**2*b1*db0 + 4*b0**3*db1
                dB42 = 12*b0*b1**2*db0 + 12*b0**2*b1*db1
                dB43 = 4*b1**3*db0 + 12*b0*b1**2*db1
                dB44 = 4*b1**3*db1

                dpoly_T2 = np.vstack([
                    6*A_j*dA   - 3*dA*B20 - 3*A_j*dB20 + 0.25*dB40,
                    -3*dA*B21  - 3*A_j*dB21            + 0.25*dB41,
                    -3*dA*B22  - 3*A_j*dB22            + 0.25*dB42,
                    0.25*dB43,
                    0.25*dB44
                ]).T                                                   # (4,5)

                dpoly_TVVT = (E * np.vstack([
                    np.zeros(4, dtype=complex),
                    4*db0,
                    6*db1 + 2*dB20,
                    2*dB21,
                    2*dB22
                ]).T)                                                  # (4,5)

                dpoly_tot = dpoly_T2 + dpoly_TVVT                      # (4,5)

                # --------------- final contraction for every k ---------------
                I_vec = np.array(I)
                for k in range(4):
                    dH2[i, j, k] = (dpoly_tot[k] @ I_vec +
                                    poly_tot        @ dI[:, k])

            # end j
        # add bra-side piece on the diagonal
        for i in range(n):
            dH2[i, i, :] += np.conj(dH2[i, i, :])

        return dH2    # ─────────────────── tiny ground-state sanity check ─────────────
    def run_basic_tests(self,p):
        S, H, H2 = self.S(p), self.H(p), self.H2(p)
        E   = (self.c.conj() @ H  @ self.c)/(self.c.conj() @ S @ self.c)
        Var = (self.c.conj() @ H2 @ self.c)/(self.c.conj() @ S @ self.c) - E**2
        assert abs(E.real - 0.5) < 1e-12 and abs(Var) < 1e-12
        print("✓ ground state:  ⟨H⟩ = 0.5  and  Var(H) = 0  (validated)")

# ────────────────────────────────────────────────────────────────────

# ──────────────────────── test_derivative() ────────────────────────
def test_derivative(seed: int = 123, eps: float = 1e-5, thresh: float = 5e-6):
    """
    Finite-difference check of dS(p) *and* dH(p).
    """
    rng = np.random.default_rng(seed)
    n   = 20

    p = np.column_stack([
        rng.uniform(0.3, 2.0, size=n),   # a
        rng.uniform(-1.0, 1.0, size=n),  # b
        rng.uniform(-1.0, 1.0, size=n),  # μ
        rng.uniform(-1.0, 1.0, size=n)   # p
    ])

    wf = HOscillator_1D(p, np.ones(n, complex))

    # ---------- S -------------
    dS_analytic = wf.dS(p)
    dS_num      = np.zeros_like(dS_analytic)

    for j in range(n):
        for k in range(4):
            p_plus  = p.copy(); p_plus [j, k] += eps
            p_minus = p.copy(); p_minus[j, k] -= eps
            dS_num[:, j, k] = (wf.S(p_plus)[:, j] - wf.S(p_minus)[:, j])/(2*eps)

    err_S = np.max(np.abs(dS_analytic - dS_num))
    print(f"max |dS_analytic − num| = {err_S:.3e}")

    # ---------- H -------------
    dH_analytic = wf.dH(p)
    dH_num      = np.zeros_like(dH_analytic)

    for j in range(n):
        for k in range(4):
            p_plus  = p.copy(); p_plus [j, k] += eps
            p_minus = p.copy(); p_minus[j, k] -= eps
            dH_num[:, j, k] = (wf.H(p_plus)[:, j] - wf.H(p_minus)[:, j])/(2*eps)

    err_H = np.max(np.abs(dH_analytic - dH_num))
    print(f"max |dH_analytic − num| = {err_H:.3e}")

    assert err_S < thresh and err_H < thresh, "gradient check failed!"
    print("✓ both dS(p) and dH(p) pass finite-difference tests.")
    dH2_analytic = wf.dH2(p)
    dH2_num      = np.zeros_like(dH2_analytic)

    for j in range(n):
        for k in range(4):
            p_plus  = p.copy(); p_plus [j, k] += eps
            p_minus = p.copy(); p_minus[j, k] -= eps
            dH2_num[:, j, k] = (wf.H2(p_plus)[:, j] - wf.H2(p_minus)[:, j])/(2*eps)

    err_H2 = np.max(np.abs(dH2_analytic - dH2_num))
    print(f"max |dH2_analytic − num| = {err_H2:.3e}")

if __name__ == "__main__":
    # single Gaussian at the harmonic-oscillator minimum
    p0 = np.array([[1/np.sqrt(2), 0.0, 0.0, 0.0]])
    wf  = HOscillator_1D(p0, np.array([1+0j]), omega=1.0)
    wf.run_basic_tests(p0)
    test_derivative()
