#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sympy as sp

I = sp.I

def coeffs(expr, x):
    expr = sp.expand(expr)
    return (sp.simplify(expr.coeff(x, 2)),
            sp.simplify(expr.coeff(x, 1)),
            sp.simplify(expr.coeff(x, 0)))

def main():
    # Symbols: all params except c_i, c_j are real
    x    = sp.symbols('x', real=True)
    ai, aj = sp.symbols('a_i a_j', real=True)
    bi, bj = sp.symbols('b_i b_j', real=True)
    mui, muj = sp.symbols('mu_i mu_j', real=True)
    pi_, pj_ = sp.symbols('p_i p_j', real=True)
    ci, cj = sp.symbols('c_i c_j')  # may be complex

    # Original Gaussians
    gi = ci*sp.exp(-(ai**2 + I*bi)*(x - mui)**2 + I*pi_*(x - mui))
    gj = cj*sp.exp(-(aj**2 + I*bj)*(x - muj)**2 + I*pj_*(x - muj))

    # Product exponent (ignore prefactors)
    E_prod = sp.expand(
        -(ai**2 + I*bi)*(x - mui)**2 + I*pi_*(x - mui)
        -(aj**2 + I*bj)*(x - muj)**2 + I*pj_*(x - muj)
    )

    # Combined parameters (general case)
    aij2 = ai**2 + aj**2
    bij  = bi + bj
    # assume aij2 != 0 for the "general case"
    muij = (ai**2*mui + aj**2*muj) / aij2
    pij  = (pi_ + pj_) + 2*(bi*(mui - muij) + bj*(muj - muij))

    # Target exponent and amplitude
    E_tgt = sp.expand(-(aij2 + I*bij)*(x - muij)**2 + I*pij*(x - muij))

    # Match x^2 and x terms
    Eprod2, Eprod1, Eprod0 = coeffs(E_prod, x)
    Etgt2,  Etgt1,  Etgt0  = coeffs(E_tgt,  x)

    assert sp.simplify(Eprod2 - Etgt2) == 0, "x^2 coefficients do not match"
    assert sp.simplify(Eprod1 - Etgt1) == 0, "x^1 coefficients do not match"

    # Constant-term difference = log(amplitude factor)
    Delta0 = sp.simplify(Eprod0 - Etgt0)

    cij = ci*cj*sp.exp(
        -(ai**2 + I*bi)*mui**2 -(aj**2 + I*bj)*muj**2
        +(aij2 + I*bij)*muij**2
        -I*(pi_*mui + pj_*muj) + I*pij*muij
    )

    # This is the decisive check: amplitude equals exp(Delta0)
    assert sp.simplify(cij - ci*cj*sp.exp(Delta0)) == 0, "Amplitude does not match"

    # Extra: verify E_prod - E_tgt really is independent of x
    assert sp.simplify(sp.diff(E_prod - E_tgt, x)) == 0, "Constant mismatch depends on x"

    # --- Degenerate edge case: a_i = a_j = 0 (pure chirps) ---
    M = sp.symbols('mu_ij_deg', real=True)  # free center
    aij2_d = 0
    bij_d  = bi + bj
    pij_d  = (pi_ + pj_) + 2*(bi*mui + bj*muj - (bi + bj)*M)

    E_prod_d = sp.expand(E_prod.subs({ai: 0, aj: 0}))
    E_tgt_d  = sp.expand(-(aij2_d + I*bij_d)*(x - M)**2 + I*pij_d*(x - M))

    Eprod2_d, Eprod1_d, Eprod0_d = coeffs(E_prod_d, x)
    Etgt2_d,  Etgt1_d,  Etgt0_d  = coeffs(E_tgt_d,  x)

    assert sp.simplify(Eprod2_d - Etgt2_d) == 0, "[degenerate] x^2 mismatch"
    assert sp.simplify(Eprod1_d - Etgt1_d) == 0, "[degenerate] x^1 mismatch"

    Delta0_d = sp.simplify(Eprod0_d - Etgt0_d)

    cij_d = ci*cj*sp.exp(
        -(I*bi)*mui**2 -(I*bj)*muj**2
        +(I*bij_d)*M**2
        -I*(pi_*mui + pj_*muj) + I*pij_d*M
    )

    assert sp.simplify(cij_d - ci*cj*sp.exp(Delta0_d)) == 0, "[degenerate] amplitude mismatch"
    assert sp.simplify(sp.diff(E_prod_d - E_tgt_d, x)) == 0, "[degenerate] constant depends on x"

    # Optional numeric smoke test (avoids heavy exact simplification of exp ratio)
    subs_numeric = {
        ai: 1.2, aj: 0.7, bi: 0.4, bj: -0.2,
        mui: -0.35, muj: 0.8, pi_: 1.1, pj_: -0.9,
        ci: 1+2*I, cj: 0.7-0.3*I, x: 0.123
    }
    lhs = sp.N((gi*gj).subs(subs_numeric))
    rhs = sp.N((cij*sp.exp(E_tgt)).subs(subs_numeric))
    assert abs(lhs - rhs) < 1e-30, "Numeric check failed in general case"

    subs_deg = {ai:0, aj:0, bi: 0.6, bj: -0.1, mui: 1.0, muj: -2.0,
                pi_: 0.5, pj_: 1.3, ci: 1-1.2*I, cj: -0.4+0.8*I, x: 2.5, M: 0.3}
    lhs_d = sp.N((gi*gj).subs(subs_deg))
    rhs_d = sp.N((cij_d*sp.exp(E_tgt_d)).subs(subs_deg))
    assert abs(lhs_d - rhs_d) < 1e-30, "Numeric check failed in degenerate case"

    print("All symbolic checks passed (general and degenerate). Numeric smoke tests passed.")

if __name__ == "__main__":
    main()
