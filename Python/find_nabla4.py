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


# ── factor out common envelope ────────────────────────────────────────────
env  = (g_i * g_j)
poly = sp.expand(product/env)
poly = sp.collect(poly, x, evaluate=False)

# build a Poly without sorting errors
poly_expr = sp.Add(*[poly[key]*key for key in poly])   # sum coeff * x^k
poly_obj = sp.Poly(poly_expr, x)
coeffs_raw = [sp.simplify(c) for c in poly_obj.all_coeffs()[::-1]]      # c₀ … c₄
print("\nRaw    FLOPs:", sum(sp.count_ops(c) for c in coeffs_raw))

# ------------------------------------------------------------------
# STEP 1 – pre-optimise each coefficient so CSE has something to bite on
# factor_terms → expand_power_exp → Horner form
# ------------------------------------------------------------------
def preshape(expr):
    expr = sp.factor_terms(expr, fraction=True)
    expr = sp.expand_power_exp(expr)
    return sp.horner(expr, x)
coeffs_pre=coeffs_raw
coeffs_pre = [preshape(c) for c in coeffs_raw]
print("After  pre-shape FLOPs:",
      sum(sp.count_ops(c) for c in coeffs_pre))

# ------------------------------------------------------------------
# STEP 2 – Common-subexpression elimination with Fortran-style temps
# ------------------------------------------------------------------
def fortran_x_symbols():
    i = 0                   # <-- start at 0 so you get x(0), x(1) …
    while True:
        yield sp.Symbol(f'in({i})')
        i += 1



def multi_cse(exprs, n_passes=1, *, optimizations=(), order='canonical'):
    """
    exprs  : sequence of SymPy expressions (kept in the final result)
    n_passes ≥ 1
    returns (all_intermediates, reduced_exprs) with len(reduced_exprs) == len(exprs)
    """
    if n_passes < 1:
        raise ValueError("n_passes must be ≥ 1")

    n_main   = len(exprs)
    sym_gen  = fortran_x_symbols()          # ONE generator
    all_ints = []

    # The working list grows with the new RHSes each round,
    # but we always keep the *last* n_main items as the main expressions.
    work = list(exprs)

    for _ in range(n_passes):
        ints, reduced = sp.cse(
            work,
            symbols       = sym_gen,
            optimizations = optimizations,
            order         = order
        )
        all_ints.extend(ints)

        # feed RHSes of new temps + current reduced back in
        work = [rhs for _, rhs in ints] + reduced
        # but only the last n_main are the 'true' outputs
        work = work[-n_main:]

    return all_ints, work                # len(work) == len(exprs)

# ──────────────────────────────────────────────────────────────────────────
# Replace your original single-pass block with the call below
# ──────────────────────────────────────────────────────────────────────────
# how many passes?  1 → original behaviour, 2 or 3 often catches more
N_PASSES =2 

cse_intermediates, cse_coeffs = multi_cse(
        coeffs_pre,
        n_passes     = N_PASSES,
        order        ='canonical'
)

# ──────────────────────────────────────────────────────────────────────────
# Pretty-print
# ──────────────────────────────────────────────────────────────────────────
print(f"After {N_PASSES}-pass CSE FLOPs:",
      sum(sp.count_ops(c) for c in cse_coeffs) +
      sum(sp.count_ops(rhs) for _, rhs in cse_intermediates))

print("\nIntermediates:")
for lhs, rhs in cse_intermediates:
    print(f"  {lhs} = {sp.simplify(rhs)}")

print("\nShortened coefficients:")
for k, coeff in enumerate(cse_coeffs):
    print(f"  c({k}) = {sp.simplify(coeff)}")          # Fortran indices









'''
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
'''