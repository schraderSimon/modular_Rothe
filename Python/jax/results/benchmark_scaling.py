import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

"""
Benchmark Rothe optimization step time as function of n (basis size) and m (potential terms).
Writes results to CSV and generates plots.
"""

import os
import time

import numpy as np

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)
# jax.config.update("jax_platform_name", "cpu")  # Force CPU

from libraries.general_wf import generalPotentialSolver
from libraries.jax_Rothe import setUpRotheErrorAndGradient_jit

# Parameters
D = 3  # dimensions
# Note: Larger n values may cause ptxas compilation issues on some GPUs
n_values = [10, 20, 40, 80, 160]  # Reduced to avoid ptxas issues
m_values = [11, 21, 31, 41, 80]
n_warmup = 3
n_runs = 10

# Output file
output_file = "benchmark_scaling_results.csv"


def build_solver(n, m, D=3):
    """Create solver with n basis functions and m Gaussian potential terms."""
    key = jax.random.PRNGKey(42)

    # Build potential string with m Gaussians
    polynomial_string = f"dimension {D}
polynomial
exponential
"
    exp_coeffs = jax.random.uniform(key, (m,), minval=0.1, maxval=2.0)
    exp_widths = jax.random.uniform(jax.random.PRNGKey(1), (m,), minval=0.1, maxval=5.0)
    centers = jax.random.uniform(jax.random.PRNGKey(2), (m, D), minval=-2.0, maxval=2.0)

    for i in range(m):
        center_str = ", ".join([f"{float(centers[i, d]):.4f}" for d in range(D)])
        polynomial_string += (
            f"{float(exp_coeffs[i])}, {float(exp_widths[i])}, [{center_str}]
"
        )

    # Random params
    params = jax.random.normal(jax.random.PRNGKey(0), (n, D, 4))
    params = params.at[:, :, 0].set(jnp.abs(params[:, :, 0]) + 0.1)

    # Coefficients
    coeffs = jax.random.normal(jax.random.PRNGKey(3), (n,), dtype=jnp.complex128)
    coeffs = coeffs / jnp.linalg.norm(coeffs)

    solver = generalPotentialSolver(params, polynomial_string)

    return solver, params, coeffs


def benchmark_optimization_step(solver, params, coeffs, n_warmup=3, n_runs=10):
    """Benchmark a single Rothe optimization step (error + gradient)."""
    SHH2 = solver.calculate_SHH2
    rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type="none")

    t = 0.0
    dt = 0.01 + 0j
    params_old = params
    params_new = params

    # Warmup
    for _ in range(n_warmup):
        _ = rothe_vg_jit(params_new, params_old, coeffs, SHH2, t, dt)

    # Benchmark
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        result = rothe_vg_jit(params_new, params_old, coeffs, SHH2, t, dt)
        # Block until ready
        result[0].block_until_ready()
        times.append(time.perf_counter() - t0)

    return np.mean(times) * 1000, np.std(times) * 1000  # ms


def main():
    print("=" * 70)
    print("Rothe Optimization Step Scaling Benchmark")
    print("=" * 70)
    print(f"D = {D}, n_warmup = {n_warmup}, n_runs = {n_runs}")
    print(f"n values: {n_values}")
    print(f"m values: {m_values}")
    print()

    # Results storage
    results = []

    # Header
    header = "n,m,time_ms,std_ms"
    print(f"{'n':>5} {'m':>5} {'time (ms)':>12} {'std (ms)':>10}")
    print("-" * 40)

    for m in m_values:
        print(f"
--- m = {m} ---")
        for n in n_values:
            try:
                solver, params, coeffs = build_solver(n, m, D)
                avg_time, std_time = benchmark_optimization_step(
                    solver, params, coeffs, n_warmup, n_runs
                )
                results.append((n, m, avg_time, std_time))
                print(f"{n:>5} {m:>5} {avg_time:>12.2f} {std_time:>10.2f}")
            except Exception as e:
                print(f"{n:>5} {m:>5} ERROR: {e}")
                results.append((n, m, float("nan"), float("nan")))

    # Write to CSV
    print(f"
Writing results to {output_file}...")
    with open(output_file, "w") as f:
        f.write(header + "
")
        for n, m, avg, std in results:
            f.write(f"{n},{m},{avg},{std}
")
    print("Done.")

    # Plot
    print("
Generating plot...")
    try:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: Time vs n for each m
        for m in m_values:
            data = [(r[0], r[2], r[3]) for r in results if r[1] == m and not np.isnan(r[2])]
            if data:
                ns, times, stds = zip(*data)
                ax1.errorbar(ns, times, yerr=stds, marker="o", label=f"m={m}", capsize=3)

        ax1.set_xlabel("n (basis functions)")
        ax1.set_ylabel("Time per optimization step (ms)")
        ax1.set_title("Rothe Optimization Step Time vs Basis Size")
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_yscale("log")

        # Plot 2: Time vs n² (should be roughly linear for matrix operations)
        for m in m_values:
            data = [(r[0], r[2], r[3]) for r in results if r[1] == m and not np.isnan(r[2])]
            if data:
                ns, times, stds = zip(*data)
                n_sq = [x**2 for x in ns]
                ax2.errorbar(n_sq, times, yerr=stds, marker="o", label=f"m={m}", capsize=3)

        ax2.set_xlabel("n² (basis functions squared)")
        ax2.set_ylabel("Time per optimization step (ms)")
        ax2.set_title("Rothe Optimization Step Time vs n²")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_file = "benchmark_scaling_plot.png"
        plt.savefig(plot_file, dpi=150)
        print(f"Plot saved to {plot_file}")
        plt.show()

    except ImportError:
        print("matplotlib not available, skipping plot")
    except Exception as e:
        print(f"Error generating plot: {e}")


if __name__ == "__main__":
    main()
