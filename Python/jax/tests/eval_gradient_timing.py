import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

"""
Benchmark test to measure which intermediate gradient is most time-consuming:
V^2, VT+TV, or T^2.

This test compares the time taken for gradient computation through each
component of the H² matrix: the kinetic squared (T²), potential squared (V²),
and the kinetic-potential cross terms (VT+TV).
"""

import time

import numpy as np

import jax
import jax.numpy as jnp
from jax import grad, jit

jax.config.update("jax_enable_x64", True)

from libraries.gaussian_potential_helpers import (
    get_exponent_subtraction_terms,
    get_gaussian_pair_terms,
)
from libraries.general_wf import dtype_complex, dtype_real, generalPotentialSolver


def calculate_Gaussian_expectation_values_autodiff(
    params, gaussian_exponential_params, linear_params
):
    """
    Same as calculate_Gaussian_expectation_values but WITHOUT custom VJP.
    Uses JAX autodiff for gradients - for comparison with analytical VJP.
    """
    n, D, _ = params.shape

    gauss_A, gauss_B, gauss_C = get_gaussian_pair_terms(params)
    exp_substr_terms = get_exponent_subtraction_terms(gauss_A, gauss_B, gauss_C)
    aK, bK, muK, pK = gaussian_exponential_params.transpose(2, 0, 1)
    alphaK = aK**2 + 1j * bK  # (m, D)

    def compute_single_contrib(carry, inputs):
        """Process one potential term at a time."""
        result = carry
        ak, muk, pk, ck = inputs

        allA = gauss_A + ak
        allB = gauss_B + 2 * (ak * muk) + 1j * pk
        allC = gauss_C - ak * muk**2 - 1j * muk * pk

        new_centers = allB / (2 * allA)
        constants = allC + new_centers**2 * allA
        sum_exponents_all = constants - 0.5 * jnp.log(allA) - 0.5 * exp_substr_terms
        sum_exponents_ij = jnp.sum(sum_exponents_all, axis=2)

        contrib = jnp.exp(sum_exponents_ij) * ck
        return result + contrib, None

    scan_inputs = (alphaK, muK, pK, linear_params)
    result, _ = jax.lax.scan(
        compute_single_contrib,
        jnp.zeros((n, n), dtype=dtype_complex),
        scan_inputs,
    )
    return result


def create_test_solver(n_gaussians: int, D: int, polynomial_string: str):
    """Create a solver with random Gaussian parameters."""
    np.random.seed(42)

    # Create random Gaussian parameters: (n, D, 4) with (a, b, mu, p)
    params = np.zeros((n_gaussians, D, 4), dtype=np.float64)
    params[:, :, 0] = np.random.uniform(0.5, 2.0, (n_gaussians, D))  # a (width)
    params[:, :, 1] = np.random.uniform(-0.5, 0.5, (n_gaussians, D))  # b (phase)
    params[:, :, 2] = np.random.uniform(-3.0, 3.0, (n_gaussians, D))  # mu (center)
    params[:, :, 3] = np.random.uniform(-1.0, 1.0, (n_gaussians, D))  # p (momentum)

    params = jnp.asarray(params, dtype=dtype_real)
    solver = generalPotentialSolver(params, polynomial_string)
    return solver, params


def benchmark_gradient_component(compute_fn, params, name, num_runs=50, warmup=10):
    """
    Benchmark the gradient of a specific component.

    Args:
        compute_fn: A JIT-compatible function that takes params and returns a matrix
        params: Initial parameters
        name: Name of the component for printing
        num_runs: Number of timing runs
        warmup: Number of warmup runs for JIT compilation

    Returns:
        Average time per gradient evaluation in milliseconds
    """

    # Create a loss function that sums all elements (to get a scalar for gradient)
    def loss(p):
        matrix = compute_fn(p)
        return jnp.sum(jnp.abs(matrix) ** 2)

    # JIT compile the gradient
    grad_fn = jit(grad(loss))

    # Warmup runs (to compile the function)
    for _ in range(warmup):
        g = grad_fn(params)
        jax.block_until_ready(g)

    # Timing runs
    start = time.perf_counter()
    for _ in range(num_runs):
        g = grad_fn(params)
        jax.block_until_ready(g)
    elapsed = time.perf_counter() - start

    avg_time_ms = (elapsed / num_runs) * 1000
    return avg_time_ms


def make_component_fn(solver, method_name, needs_t=False):
    """Create a function that computes a specific matrix component."""

    def compute(p):
        solver.update_parameters(p)
        solver.S = None  # Force S recalculation
        if needs_t:
            return getattr(solver, method_name)(t=0)
        else:
            return getattr(solver, method_name)()

    return compute


def benchmark_vjp_speedup(
    n_gaussians: int = 20,
    D: int = 3,
    n_potential_gaussians: int = 10,
    num_runs: int = 50,
    warmup: int = 10,
):
    """
    Compare gradient speed with analytical VJP vs autodiff for Gaussian potentials.
    """
    from libraries.gaussian_potential_helpers import (
        calculate_Gaussian_expectation_values,
        calculate_Gaussian_expectation_values_batched,
        calculate_Gaussian_expectation_values_fast,
        compute_batch_size,
    )

    print(f"\n{'='*70}")
    print(f"VJP Comparison: Original vs Fast Vectorized vs Batched vs Autodiff")
    print(
        f"n={n_gaussians} basis Gaussians, D={D} dims, m={n_potential_gaussians} potential Gaussians"
    )
    print(f"{'='*70}")

    # Create random params
    np.random.seed(42)
    params = np.zeros((n_gaussians, D, 4), dtype=np.float64)
    params[:, :, 0] = np.random.uniform(0.5, 2.0, (n_gaussians, D))
    params[:, :, 1] = np.random.uniform(-0.5, 0.5, (n_gaussians, D))
    params[:, :, 2] = np.random.uniform(-3.0, 3.0, (n_gaussians, D))
    params[:, :, 3] = np.random.uniform(-1.0, 1.0, (n_gaussians, D))
    params = jnp.asarray(params, dtype=dtype_real)

    # Create potential params (fixed, not differentiated)
    widths = np.logspace(-1, 2, n_potential_gaussians)
    pot_params = np.zeros((n_potential_gaussians, D, 4), dtype=np.float64)
    for i, w in enumerate(widths):
        pot_params[i, :, 0] = w
    pot_params = jnp.asarray(pot_params, dtype=dtype_real)
    lin_params = jnp.array(
        [(-1) ** i / (i + 1) for i in range(n_potential_gaussians)], dtype=dtype_complex
    )

    # Test V (single Gaussian potential layer)
    print("\\n--- V (Gaussian potential) ---")

    def loss_vjp(p):
        result = calculate_Gaussian_expectation_values(p, pot_params, lin_params)
        return jnp.sum(jnp.abs(result) ** 2)

    def loss_vjp_fast(p):
        result = calculate_Gaussian_expectation_values_fast(p, pot_params, lin_params)
        return jnp.sum(jnp.abs(result) ** 2)

    def loss_vjp_batched(p):
        result = calculate_Gaussian_expectation_values_batched(p, pot_params, lin_params, 500)
        return jnp.sum(jnp.abs(result) ** 2)

    def loss_autodiff(p):
        result = calculate_Gaussian_expectation_values_autodiff(p, pot_params, lin_params)
        return jnp.sum(jnp.abs(result) ** 2)

    grad_vjp = jit(grad(loss_vjp))
    grad_vjp_fast = jit(grad(loss_vjp_fast))
    grad_vjp_batched = jit(grad(loss_vjp_batched))
    grad_autodiff = jit(grad(loss_autodiff))

    # Verify correctness first
    g1 = grad_vjp(params)
    g2 = grad_vjp_fast(params)
    g3 = grad_vjp_batched(params)
    g4 = grad_autodiff(params)

    err_fast_vs_orig = jnp.max(jnp.abs(g1 - g2))
    err_batched_vs_orig = jnp.max(jnp.abs(g1 - g3))
    err_fast_vs_auto = jnp.max(jnp.abs(g2 - g4))
    print(f"  Max error (fast vs original VJP):    {err_fast_vs_orig:.2e}")
    print(f"  Max error (batched vs original VJP): {err_batched_vs_orig:.2e}")
    print(f"  Max error (fast vs autodiff):        {err_fast_vs_auto:.2e}")

    # Warmup
    for _ in range(warmup):
        jax.block_until_ready(grad_vjp(params))
        jax.block_until_ready(grad_vjp_fast(params))
        jax.block_until_ready(grad_vjp_batched(params))
        jax.block_until_ready(grad_autodiff(params))

    # Time original VJP
    start = time.perf_counter()
    for _ in range(num_runs):
        jax.block_until_ready(grad_vjp(params))
    time_vjp = (time.perf_counter() - start) / num_runs * 1000

    # Time fast VJP
    start = time.perf_counter()
    for _ in range(num_runs):
        jax.block_until_ready(grad_vjp_fast(params))
    time_vjp_fast = (time.perf_counter() - start) / num_runs * 1000

    # Time batched VJP
    start = time.perf_counter()
    for _ in range(num_runs):
        jax.block_until_ready(grad_vjp_batched(params))
    time_vjp_batched = (time.perf_counter() - start) / num_runs * 1000

    # Time autodiff
    start = time.perf_counter()
    for _ in range(num_runs):
        jax.block_until_ready(grad_autodiff(params))
    time_autodiff = (time.perf_counter() - start) / num_runs * 1000

    print(f"  Original VJP:    {time_vjp:8.3f} ms")
    print(
        f"  Fast VJP:        {time_vjp_fast:8.3f} ms  (speedup vs original: {time_vjp/time_vjp_fast:.2f}x)"
    )
    print(
        f"  Batched VJP:     {time_vjp_batched:8.3f} ms  (speedup vs original: {time_vjp/time_vjp_batched:.2f}x)"
    )
    print(
        f"  Autodiff:        {time_autodiff:8.3f} ms  (speedup vs autodiff: {time_autodiff/time_vjp_fast:.2f}x)"
    )

    # Test V² (squared Gaussian potential - more terms)
    print("\\n--- V² (squared Gaussian potential) ---")

    # Create squared potential params
    polynomial_string = create_gaussian_potential_string(D, n_potential_gaussians)
    solver = generalPotentialSolver(params, polynomial_string)

    sq_params = solver.exponential_squared_params
    sq_lin = solver.linear_exponential_squared
    n_sq = sq_params.shape[0] if sq_params is not None else 0
    print(f"  (V² uses {n_sq} squared Gaussian terms)")

    # Estimate memory for fast VJP
    if n_sq > 0:
        bytes_per_tensor = n_sq * n_gaussians * n_gaussians * D * 16
        est_mem_gb = bytes_per_tensor * 15 / 1e9
        batch_size = compute_batch_size(n_sq, n_gaussians, D, 500_000_000)  # 500MB limit
        print(
            f"  (Est. fast VJP memory: {est_mem_gb:.2f} GB, batched will use batch_size={batch_size})"
        )

    if sq_params is not None and n_sq > 0:

        def loss_vjp_sq(p):
            result = calculate_Gaussian_expectation_values(p, sq_params, sq_lin)
            return jnp.sum(jnp.abs(result) ** 2)

        def loss_vjp_fast_sq(p):
            result = calculate_Gaussian_expectation_values_fast(p, sq_params, sq_lin)
            return jnp.sum(jnp.abs(result) ** 2)

        def loss_vjp_batched_sq(p):
            result = calculate_Gaussian_expectation_values_batched(p, sq_params, sq_lin, 500)
            return jnp.sum(jnp.abs(result) ** 2)

        def loss_autodiff_sq(p):
            result = calculate_Gaussian_expectation_values_autodiff(p, sq_params, sq_lin)
            return jnp.sum(jnp.abs(result) ** 2)

        grad_vjp_sq = jit(grad(loss_vjp_sq))
        grad_vjp_fast_sq = jit(grad(loss_vjp_fast_sq))
        grad_vjp_batched_sq = jit(grad(loss_vjp_batched_sq))
        grad_autodiff_sq = jit(grad(loss_autodiff_sq))

        # Verify correctness
        g1 = grad_vjp_sq(params)
        g2 = grad_vjp_fast_sq(params)
        g3 = grad_vjp_batched_sq(params)
        g4 = grad_autodiff_sq(params)

        err_fast_vs_orig = jnp.max(jnp.abs(g1 - g2))
        err_batched_vs_orig = jnp.max(jnp.abs(g1 - g3))
        err_fast_vs_auto = jnp.max(jnp.abs(g2 - g4))
        print(f"  Max error (fast vs original VJP):    {err_fast_vs_orig:.2e}")
        print(f"  Max error (batched vs original VJP): {err_batched_vs_orig:.2e}")
        print(f"  Max error (fast vs autodiff):        {err_fast_vs_auto:.2e}")

        # Warmup
        for _ in range(warmup):
            jax.block_until_ready(grad_vjp_sq(params))
            jax.block_until_ready(grad_vjp_fast_sq(params))
            jax.block_until_ready(grad_vjp_batched_sq(params))
            jax.block_until_ready(grad_autodiff_sq(params))

        # Time original VJP
        start = time.perf_counter()
        for _ in range(num_runs):
            jax.block_until_ready(grad_vjp_sq(params))
        time_vjp_sq = (time.perf_counter() - start) / num_runs * 1000

        # Time fast VJP
        start = time.perf_counter()
        for _ in range(num_runs):
            jax.block_until_ready(grad_vjp_fast_sq(params))
        time_vjp_fast_sq = (time.perf_counter() - start) / num_runs * 1000

        # Time batched VJP
        start = time.perf_counter()
        for _ in range(num_runs):
            jax.block_until_ready(grad_vjp_batched_sq(params))
        time_vjp_batched_sq = (time.perf_counter() - start) / num_runs * 1000

        # Time autodiff
        start = time.perf_counter()
        for _ in range(num_runs):
            jax.block_until_ready(grad_autodiff_sq(params))
        time_autodiff_sq = (time.perf_counter() - start) / num_runs * 1000

        print(f"  Original VJP:    {time_vjp_sq:8.3f} ms")
        print(
            f"  Fast VJP:        {time_vjp_fast_sq:8.3f} ms  (speedup vs original: {time_vjp_sq/time_vjp_fast_sq:.2f}x)"
        )
        print(
            f"  Batched VJP:     {time_vjp_batched_sq:8.3f} ms  (speedup vs original: {time_vjp_sq/time_vjp_batched_sq:.2f}x)"
        )
        print(
            f"  Autodiff:        {time_autodiff_sq:8.3f} ms  (speedup vs autodiff: {time_autodiff_sq/time_vjp_fast_sq:.2f}x)"
        )

    print(f"\\n{'='*70}")


def create_gaussian_potential_string(D: int, n_potential_gaussians: int = 10):
    """
    Create a potential string with Gaussian (exponential) terms.
    This mimics a Coulomb-like potential approximated by Gaussians.
    """
    # Use a range of widths similar to the Hydrogen potential approximation
    widths = np.logspace(-1, 2, n_potential_gaussians)

    lines = [f"dimension {D}", "exponential"]
    center = ", ".join(["0.0"] * D)

    for i, w in enumerate(widths):
        # Alternate signs to create an interesting potential
        coeff = (-1) ** i * (1.0 / (i + 1))
        lines.append(f"{coeff}, {w}, [{center}]")

    return "\n".join(lines)


def create_polynomial_potential_string(D: int):
    """Create a polynomial potential string (Henon-Heiles like for 2D/3D)."""
    if D == 2:
        d = """
           dimension 2
           polynomial
           x0x0: 0.5
           x1x1: 0.5
           x0x0x1: 0.111803
           x0x1x1: -0.0372678
            """
        return d
    elif D == 3:
        d = """
    dimension 3
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x2x2: 0.5
    x0x0x1: 0.1
    x1x1x2: 0.1
    """
    else:
        # Generic harmonic oscillator for other dimensions
        lines = [f"dimension {D}", "polynomial"]
        for d in range(D):
            lines.append(f"x{d}x{d}: 0.5")
        return "\n".join(lines)


def create_mixed_potential_string(D: int, n_potential_gaussians: int = 5):
    """Create a potential with both polynomial and Gaussian terms."""
    if D == 2:
        lines = [
            "dimension 2",
            "polynomial",
            "x0x0: 0.5",
            "x1x1: 0.5",
            "exponential",
        ]
    elif D == 3:
        lines = [
            "dimension 3",
            "polynomial",
            "x0x0: 0.5",
            "x1x1: 0.5",
            "x2x2: 0.5",
            "exponential",
        ]
    else:
        lines = [f"dimension {D}", "polynomial"]
        for d in range(D):
            lines.append(f"x{d}x{d}: 0.5")
        lines.append("exponential")

    widths = np.logspace(-1, 2, n_potential_gaussians)
    center = ", ".join(["0.0"] * D)
    for i, w in enumerate(widths):
        coeff = (-1) ** i * (1.0 / (i + 1))
        lines.append(f"{coeff}, {w}, [{center}]")

    return "\n".join(lines)


def run_benchmark(
    n_gaussians: int = 20,
    D: int = 3,
    num_runs: int = 50,
    warmup: int = 10,
    potential_type: str = "gaussian",
    n_potential_gaussians: int = 10,
):
    """
    Run the full benchmark comparing V², VT+TV, and T² gradient times.

    Args:
        n_gaussians: Number of basis Gaussians
        D: Number of dimensions
        num_runs: Number of timing runs
        warmup: Number of warmup runs
        potential_type: "gaussian", "polynomial", or "mixed"
        n_potential_gaussians: Number of Gaussians in the potential (for gaussian/mixed types)
    """
    print("=" * 70)
    print(f"Gradient Timing Benchmark: n={n_gaussians} Gaussians, D={D} dimensions")
    print(f"Potential type: {potential_type}", end="")
    if potential_type in ("gaussian", "mixed"):
        print(f" ({n_potential_gaussians} potential Gaussians)")
    else:
        print()
    print(f"{'='*70}")

    # Create potential string based on type
    if potential_type == "gaussian":
        polynomial_string = create_gaussian_potential_string(D, n_potential_gaussians)
    elif potential_type == "polynomial":
        polynomial_string = create_polynomial_potential_string(D)
    elif potential_type == "mixed":
        polynomial_string = create_mixed_potential_string(D, n_potential_gaussians)
    else:
        raise ValueError(f"Unknown potential_type: {potential_type}")

    # Create initial solver and params
    _, params = create_test_solver(n_gaussians, D, polynomial_string)

    # Create a reference solver ONCE outside of JIT
    # This parses the potential string and sets up exponential_params, etc.
    solver = generalPotentialSolver(params, polynomial_string)

    # Create compute functions that update params on the existing solver
    # This avoids re-parsing the potential string inside JIT
    compute_Tsq = make_component_fn(solver, "calculate_Tsq", needs_t=False)
    compute_V2 = make_component_fn(solver, "calculate_V2", needs_t=True)
    compute_TVVT = make_component_fn(solver, "calculate_TVVT", needs_t=True)
    compute_H2 = make_component_fn(solver, "calculate_H2", needs_t=True)
    compute_S = make_component_fn(solver, "calculate_S", needs_t=False)
    compute_H = make_component_fn(solver, "calculate_H", needs_t=True)

    print(f"Running {num_runs} iterations with {warmup} warmup runs each...")

    # Benchmark T² (kinetic squared)
    time_Tsq = benchmark_gradient_component(compute_Tsq, params, "T²", num_runs, warmup)
    print(f"  T² (kinetic squared):     {time_Tsq:8.3f} ms")

    # Benchmark V² (potential squared)
    time_V2 = benchmark_gradient_component(compute_V2, params, "V²", num_runs, warmup)
    print(f"  V² (potential squared):   {time_V2:8.3f} ms")

    # Benchmark VT+TV (kinetic-potential cross terms)
    time_TVVT = benchmark_gradient_component(compute_TVVT, params, "VT+TV", num_runs, warmup)
    print(f"  VT+TV (cross terms):      {time_TVVT:8.3f} ms")

    # Benchmark full H² for comparison
    time_H2 = benchmark_gradient_component(compute_H2, params, "H²", num_runs, warmup)
    print(f"  H² (full):                {time_H2:8.3f} ms")

    # Also benchmark S and H for reference
    time_S = benchmark_gradient_component(compute_S, params, "S", num_runs, warmup)
    print(f"  S (overlap):              {time_S:8.3f} ms")

    time_H = benchmark_gradient_component(compute_H, params, "H", num_runs, warmup)
    print(f"  H (Hamiltonian):          {time_H:8.3f} ms")

    # Summary
    total_components = time_Tsq + time_V2 + time_TVVT
    print()
    print(f"  Sum of components:        {total_components:8.3f} ms")
    print()

    # Find the most time-consuming component
    times = {
        "T² (kinetic squared)": time_Tsq,
        "V² (potential squared)": time_V2,
        "VT+TV (cross terms)": time_TVVT,
    }
    slowest = max(times, key=times.get)

    print(f"{'='*70}")
    print(f"RESULT: The most time-consuming gradient is: {slowest}")
    print(
        f"        ({times[slowest]:.3f} ms, "
        f"{times[slowest]/total_components*100:.1f}% of H² components)"
    )
    print(f"{'='*70}")

    # Relative timings
    print("Relative timings (normalized to T²):")

    print(f"  T²:    {time_Tsq/time_Tsq:6.2f}x")
    print(f"  V²:    {time_V2/time_Tsq:6.2f}x")
    print(f"  VT+TV: {time_TVVT/time_Tsq:6.2f}x")

    return {
        "T²": time_Tsq,
        "V²": time_V2,
        "VT+TV": time_TVVT,
        "H²": time_H2,
        "S": time_S,
        "H": time_H,
    }


def run_scaling_study(
    dimensions=(2, 3),
    n_range=(5, 10, 20, 30),
    num_runs=30,
    warmup=5,
    potential_type="gaussian",
    n_potential_gaussians=10,
):
    """
    Run benchmarks across different problem sizes to understand scaling.
    """
    print("" + "=" * 70)

    print(f"SCALING STUDY: Gradient timing vs problem size")
    print(f"Potential type: {potential_type}", end="")
    if potential_type in ("gaussian", "mixed"):
        print(f" ({n_potential_gaussians} potential Gaussians)")
    else:
        print()
    print("=" * 70)

    results = []

    for D in dimensions:
        print(f"--- Dimension D = {D} ---")
        for n in n_range:
            print(f"n = {n} Gaussians:")
            times = run_benchmark(
                n_gaussians=n,
                D=D,
                num_runs=num_runs,
                warmup=warmup,
                potential_type=potential_type,
                n_potential_gaussians=n_potential_gaussians,
            )
            results.append({"n": n, "D": D, **times})

    # Print summary table
    print("" + "=" * 70)

    print("SUMMARY TABLE (times in ms)")
    print("=" * 70)
    print(f"{'D':>3} {'n':>5} {'T²':>10} {'V²':>10} {'VT+TV':>10} {'H²':>10} {'Slowest':>12}")
    print("-" * 70)

    for r in results:
        slowest = (
            "VT+TV"
            if r["VT+TV"] > max(r["T²"], r["V²"])
            else ("V²" if r["V²"] > r["T²"] else "T²")
        )
        print(
            f"{r['D']:>3} {r['n']:>5} {r['T²']:>10.3f} {r['V²']:>10.3f} "
            f"{r['VT+TV']:>10.3f} {r['H²']:>10.3f} {slowest:>12}"
        )

    return results


if __name__ == "__main__":
    import sys

    # Parse command line arguments
    args = sys.argv[1:]
    potential_type = "gaussian"  # Default to Gaussian potential
    n_potential_gaussians = 10

    # Check for potential type flags
    if "--polynomial" in args:
        potential_type = "polynomial"
        args.remove("--polynomial")
    elif "--mixed" in args:
        potential_type = "mixed"
        args.remove("--mixed")
    elif "--gaussian" in args:
        potential_type = "gaussian"
        args.remove("--gaussian")

    # Check for number of potential Gaussians
    for i, arg in enumerate(args):
        if arg.startswith("--n-pot="):
            n_potential_gaussians = int(arg.split("=")[1])
            args.remove(arg)
            break

    if "--scaling" in args:
        # Run full scaling study
        run_scaling_study(
            potential_type=potential_type, n_potential_gaussians=n_potential_gaussians
        )
    elif "--compare" in args:
        # Compare all potential types
        print("" + "=" * 70)

        print("COMPARISON OF POTENTIAL TYPES")
        print("=" * 70)

        for ptype in ["polynomial", "gaussian", "mixed"]:
            run_benchmark(
                n_gaussians=20,
                D=3,
                num_runs=50,
                warmup=10,
                potential_type=ptype,
                n_potential_gaussians=10,
            )
    elif "--vjp" in args:
        # Compare analytical VJP vs autodiff
        benchmark_vjp_speedup(
            n_gaussians=20,
            D=3,
            n_potential_gaussians=n_potential_gaussians,
            num_runs=50,
            warmup=10,
        )
    else:
        # Run single benchmark with default parameters
        run_benchmark(
            n_gaussians=20,
            D=3,
            num_runs=50,
            warmup=10,
            potential_type=potential_type,
            n_potential_gaussians=n_potential_gaussians,
        )

        print("" + "-" * 70)

        print("Options:")
        print("  --scaling      Run full scaling study across problem sizes")
        print("  --compare      Compare polynomial, gaussian, and mixed potentials")
        print("  --vjp          Compare analytical VJP vs autodiff for Gaussian V/V²")
        print("  --polynomial   Use polynomial potential only")
        print("  --gaussian     Use Gaussian potential only (default)")
        print("  --mixed        Use both polynomial and Gaussian potential")
        print("  --n-pot=N      Number of Gaussians in the potential (default: 10)")
        print("-" * 70)
