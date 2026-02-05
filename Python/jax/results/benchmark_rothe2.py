import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

"""
Comprehensive benchmark of Rothe optimization components.
Profiles: S, T, V, T², V², TV+VT, and their gradients.

Uses the actual Rothe error computation workflow which is properly differentiable.
"""

import time

import numpy as np

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from libraries.general_wf import generalPotentialSolver
from libraries.jax_Rothe import setUpRotheErrorAndGradient_jit

# Problem size
n = 40  # basis functions
D = 3  # dimensions
m = 20  # potential terms

print(f"Problem size: n={n}, D={D}, m={m}")
print()

# Build potential string
polynomial_string = "dimension 3
polynomial
exponential
"
key = jax.random.PRNGKey(42)
exp_coeffs = jax.random.uniform(key, (m,), minval=0.1, maxval=2.0)
exp_widths = jax.random.uniform(jax.random.PRNGKey(1), (m,), minval=0.1, maxval=5.0)
for i in range(m):
    polynomial_string += f"{float(exp_coeffs[i])}, {float(exp_widths[i])}, [0.0, 0.0, 0.0]
"

# Random params
params = jax.random.normal(jax.random.PRNGKey(0), (n, D, 4))
params = params.at[:, :, 0].set(jnp.abs(params[:, :, 0]) + 0.1)

# Coefficients
coeffs = jax.random.normal(jax.random.PRNGKey(3), (n,), dtype=jnp.complex128)
coeffs = coeffs / jnp.linalg.norm(coeffs)

print("Creating solver...")
solver = generalPotentialSolver(params, polynomial_string)
SHH2 = solver.calculate_SHH2
print("Done.")
print()


def benchmark(name, fn, n_warmup=3, n_runs=10, return_result=False):
    for _ in range(n_warmup):
        result = fn()
        if hasattr(result, "block_until_ready"):
            result.block_until_ready()
        elif isinstance(result, tuple):
            for r in result:
                if hasattr(r, "block_until_ready"):
                    r.block_until_ready()

    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        result = fn()
        if hasattr(result, "block_until_ready"):
            result.block_until_ready()
        elif isinstance(result, tuple):
            for r in result:
                if hasattr(r, "block_until_ready"):
                    r.block_until_ready()
        times.append(time.perf_counter() - t0)

    avg = 1000 * np.mean(times)
    std = 1000 * np.std(times)
    print(f"{name:25s}: {avg:8.2f} ± {std:.2f} ms")
    if return_result:
        return avg, result
    return avg


print("=" * 60)
print("FORWARD PASS (matrix computation)")
print("=" * 60)

# Benchmark intermediates breakdown
results = {}


def setup_basic_intermediates():
    """Only Gamma, tilde_k, tilde_y, etc. - NOT moments"""
    solver.S = None  # Reset
    alpha = solver.a**2 + 1j * solver.b
    ai, aj = alpha.conj()[:, None, :], alpha[None, :, :]
    Gamma = ai + aj
    pi, pj = solver.p[:, None, :], solver.p[None, :, :]
    mui = solver.mu[:, None, :] - 1j * pi / (2 * ai)
    muj = solver.mu[None, :, :] + 1j * pj / (2 * aj)
    exponent_contrib = -(pi.conj() ** 2) / (4 * ai) - (pj**2) / (4 * aj)
    tilde_k = (ai * aj) / Gamma
    tilde_y = mui - muj
    tysq = tilde_y**2
    tksq = tilde_k**2
    return Gamma, tilde_k, tilde_y, tysq, tksq


results["basic_intermediates"] = benchmark("basic intermediates", setup_basic_intermediates)

# Moments only needed for polynomial potentials
if solver.polynomial is not None:

    def setup_moments_only():
        """Just the moments computation"""
        return solver.calculate_all_moments()

    results["moments"] = benchmark("moments (polynomial)", setup_moments_only)
else:
    results["moments"] = 0.0
    print(f"{'moments (polynomial)':25s}:     0.00 ± 0.00 ms (no polynomial)")


def setup_and_calculate_S():
    solver.setUpIntermediates()
    return solver.calculate_S()


results["intermediates+S"] = benchmark("full intermediates + S", setup_and_calculate_S)
results["S"] = benchmark("S (overlap)", lambda: solver.calculate_S())

results["T"] = benchmark("T (kinetic)", lambda: solver.calculate_T())
results["T²"] = benchmark("T² (kinetic²)", lambda: solver.calculate_Tsq())
results["V"] = benchmark("V (potential)", lambda: solver.calculate_V())
results["V²"] = benchmark("V² (potential²)", lambda: solver.calculate_V2())
results["TV+VT"] = benchmark("TV+VT (cross)", lambda: solver.calculate_TVVT())
results["H"] = benchmark("H (full)", lambda: solver.calculate_H())
results["H²"] = benchmark("H² (full)", lambda: solver.calculate_H2())
results["SHH2"] = benchmark("SHH2 (all)", lambda: solver.calculate_SHH2())


# Benchmark full update_parameters + SHH2 (what actually happens in optimization)
def full_update_and_SHH2():
    solver.update_parameters(params)
    return solver.calculate_SHH2()


results["update+SHH2"] = benchmark("update_params + SHH2", full_update_and_SHH2)

print()
print("=" * 60)
print("ROTHE ERROR + GRADIENT (actual optimization step)")
print("=" * 60)

# Use the actual Rothe error and gradient computation
rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type="none")

t = 0.0
dt = 0.01 + 0j  # complex dt for imaginary time
params_old = params
params_new = params  # In practice these differ, but for benchmarking same shape

# Warmup
_ = rothe_vg_jit(params_new, params_old, coeffs, SHH2, t, dt)

results["Rothe_vg"] = benchmark(
    "Rothe error + grad", lambda: rothe_vg_jit(params_new, params_old, coeffs, SHH2, t, dt)
)

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)

fwd = results["SHH2"]
rothe_total = results["Rothe_vg"]

print()
print(f"Single SHH2 forward:      {fwd:8.2f} ms")
print(f"With intermediates:       {results['update+SHH2']:8.2f} ms")
print(f"Rothe error + gradient:   {rothe_total:8.2f} ms")
print()
print("Intermediates breakdown:")
print(f"  basic (k,y,Gamma):   {results['basic_intermediates']:6.2f} ms")
print(f"  moments (poly only): {results['moments']:6.2f} ms")
print()
print("Forward pass breakdown (uses precomputed intermediates):")
for key in ["S", "T", "T²", "V", "V²", "TV+VT"]:
    pct = 100 * results[key] / fwd if fwd > 0 else 0
    print(f"  {key:12s}: {results[key]:6.2f} ms ({pct:5.1f}%)")

print()
print("Note: Rothe uses 2n×2n matrices (params_old + params_new concatenated)")
print(f"So expect ~4x the SHH2 time = {4*fwd:.0f} ms, plus gradient overhead")
