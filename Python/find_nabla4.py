#!/usr/bin/env python3
"""
(∇² g_i)(∇² g_j)  factorised as  env(x) * P(x)
with   P(x) = Σ_{k=0}^{4} c_k x^k
"""

import sympy as sp
sp.init_printing(use_unicode=True)

# ── symbols ────────────────────────────────────────────────────────────────
x      = sp.symbols('x', real=True)

cA_i, A_j = sp.symbols('bA_i_c bA_j', complex=True)
mu_i, mu_j = sp.symbols('mu_i mu_j', real=True)
p_i,  p_j  = sp.symbols('p_i  p_j',  real=True)
N_i, N_j,omega   = sp.symbols('N_i N_j,omega',   positive=True, real=True)

# ── Gaussians ──────────────────────────────────────────────────────────────
g_i = N_i * sp.exp(-cA_i*(x - mu_i)**2 - sp.I*p_i*(x - mu_i))
g_j = N_j * sp.exp( -A_j*(x - mu_j)**2 + sp.I*p_j*(x - mu_j))

# ── Laplacians ─────────────────────────────────────────────────────────────
lap_gi = sp.diff(g_i, x, 2)
lap_gj = sp.diff(g_j, x, 2)
product = sp.simplify(lap_gi * lap_gj)

# ── introduce α and β before expanding ────────────────────────────────────
alpha, beta = sp.symbols('alpha beta')
subs = {
    alpha: cA_i + A_j,
    beta : 2*(cA_i*mu_i + A_j*mu_j) + sp.I*(p_j - p_i)
}
product = product.subs(subs)

# ── factor out common envelope ────────────────────────────────────────────
env  = (g_i * g_j).subs(subs)          # same envelope, nicer parameters
poly = sp.expand(product/env)
poly = sp.collect(poly, x, evaluate=False)

# build a Poly without sorting errors
poly_expr = sp.Add(*[poly[key]*key for key in poly])   # sum coeff * x^k
poly_obj  = sp.Poly(poly_expr, x)

# ── output ────────────────────────────────────────────────────────────────
print("\nEnvelope  e(x) =")
sp.pprint(env, wrap_line=False)

print("\nPolynomial coefficients:")
for k, coeff in enumerate(poly_obj.all_coeffs()[::-1]):   # c_0 … c_4
    print(f"  c_{k} =", sp.simplify(coeff))

print("\nLaTeX  P(x) =")
print(sp.latex(poly_obj.as_expr()))


"""Kinetic times potential"""
product= sp.simplify(lap_gi * g_j * x**2)
poly=sp.expand(product/env)
poly = sp.collect(poly, x, evaluate=False)
poly_expr = sp.Add(*[poly[key]*key for key in poly])   # sum coeff * x^k
poly_obj  = sp.Poly(poly_expr, x)
print("\n Kinetic times potential coefficients:")
for k, coeff in enumerate(poly_obj.all_coeffs()[::-1]):   # c_0 … c_4
    print(f"  c{k} =", sp.simplify(coeff))

"""Potential times kinetic"""
product= sp.simplify(g_i * lap_gj * x**2)
poly=sp.expand(product/env)
poly = sp.collect(poly, x, evaluate=False)
poly_expr = sp.Add(*[poly[key]*key for key in poly])   # sum coeff * x^k
poly_obj  = sp.Poly(poly_expr, x)
print("\n Kinetic times potential coefficients:")
for k, coeff in enumerate(poly_obj.all_coeffs()[::-1]):   # c_0 … c_4
    print(f"  c{k} =", sp.simplify(coeff))