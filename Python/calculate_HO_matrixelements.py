#!/usr/bin/env python3
# check_overlap_grads.py
#
# Overlap S_ij of 1-D complex Gaussians and analytic gradients.
# Confirms that compact “hand” formulas equal SymPy’s raw derivatives.

import sympy as sp
import random, math

# ────────────────────────── 1. symbols ──────────────────────────
a_i, a_j = sp.symbols('a_i a_j', positive=True, real=True)
b_i, b_j = sp.symbols('b_i b_j', real=True)
mu_i, mu_j = sp.symbols('mu_i mu_j', real=True)
p_i,  p_j  = sp.symbols('p_i p_j',  real=True)

# ────────────────────────── 2. combos ───────────────────────────
A_i, A_j = a_i**2 + sp.I*b_i, a_j**2 + sp.I*b_j

alpha = sp.conjugate(A_i) + A_j
beta  = 2*(sp.conjugate(A_i)*mu_i + A_j*mu_j) + sp.I*(p_j - p_i)
gamma = -(sp.conjugate(A_i)*mu_i**2 + A_j*mu_j**2) \
        + sp.I*(p_i*mu_i - p_j*mu_j)

N_i = (2*a_i**2/sp.pi)**sp.Rational(1, 4)
N_j = (2*a_j**2/sp.pi)**sp.Rational(1, 4)

S_ij = N_i * N_j * sp.sqrt(sp.pi/alpha) * sp.exp(gamma + beta**2/(4*alpha))

# ───────────────────────── 3. raw ∂ lnS / ∂p  ───────────────────
ket_params = (a_j, b_j, mu_j, p_j)
dlogS_code = {p: sp.diff(sp.log(S_ij), p) for p in ket_params}

# ─────────────────────── 4. compact analytic forms ──────────────
expr_manual = {
    a_j: 1/(2*a_j) - a_j/alpha + 2*a_j*mu_j*beta/alpha
          - 2*a_j*mu_j**2 - a_j*beta**2/(2*alpha**2),

    b_j: sp.I*(-mu_j**2 - 1/(2*alpha) + mu_j*beta/alpha
               - beta**2/(4*alpha**2)),

    mu_j: A_j*(beta/alpha - 2*mu_j) - sp.I*p_j,

    p_j: sp.I*(beta/(2*alpha) - mu_j)
}

# ─────────────────── 5. algebraic identity check ────────────────
print("\nAlgebraic check  (should all be 0):\n")
for p in ket_params:
    diff = sp.simplify(dlogS_code[p] - expr_manual[p])
    print(f"∂lnS/∂{p}  →", diff)

# ───────────────────── 6. numerical spot-check ──────────────────
rng = random.Random(42)
subs_num = {
    a_i: 0.8,
    a_j: 1.3,
    b_i: -0.4,
    b_j: 0.55,
    mu_i: 0.15,
    mu_j: -0.6,
    p_i: 0.25,
    p_j: -0.17
}

max_err = 0.0
for p in ket_params:
    val_code   = complex(sp.N(dlogS_code[p].subs(subs_num)))
    val_manual = complex(sp.N(expr_manual[p].subs(subs_num)))
    err = abs(val_code - val_manual)
    max_err = max(max_err, err)
    print(f"{p}: code = {val_code: .6g}   manual = {val_manual: .6g}   diff = {err:.3e}")

print("\nmaximum abs error across the four params:", max_err)
assert max_err < 1e-12, "derivatives disagree!"
