import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

"""
Timing comparison of gradient implementations for Gaussian potential.

Tests:
1. Naive autodiff (no custom_vjp)
2. Fast (fully vectorized custom VJP)
3. Batched (10 batches, regardless of memory)
"""

import time

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from libraries.gaussian_potential_vjp_fast import (
    _calculate_Gaussian_expectation_values_raw,
    calculate_Gaussian_expectation_values_batched,
    calculate_Gaussian_expectation_values_fast,
)


def time_gradient(grad_fn, params, num_runs=5):
    """Time a gradient function, return mean time after warmup."""
    # Warmup / compile
    grad = grad_fn(params)
    grad.block_until_ready()

    times = []
    for _ in range(num_runs):
        t0 = time.perf_counter()
        grad = grad_fn(params)
        grad.block_until_ready()
        times.append(time.perf_counter() - t0)

    return sum(times) / len(times), grad


def run_benchmark(n, m, D=3):
    """Run benchmark for a given problem size."""
    print(f"\n{'='*70}")
    print(f"Problem size: n={n}, m={m}, D={D}")
    print(f"{'='*70}")

    # Create test data
    key = jax.random.PRNGKey(42)
    key1, key2, key3 = jax.random.split(key, 3)

    params = jax.random.normal(key1, (n, D, 4))
    params = params.at[:, :, 0].set(jnp.abs(params[:, :, 0]) + 0.5)

    pot_params = jax.random.normal(key2, (m, D, 4))
    pot_params = pot_params.at[:, :, 0].set(jnp.abs(pot_params[:, :, 0]) + 0.5)

    linear_params = jax.random.normal(key3, (m,)) + 0j

    # Calculate memory limit for exactly 10 batches
    # Memory per batch ≈ 15 * batch_size * n² * D * 32 bytes (complex128)
    # batch_size = m / 10, so memory_limit = 15 * (m/10) * n² * D * 32 / 1e6 MB
    batch_size = m // 10
    if batch_size < 1:
        batch_size = 1
    memory_limit_mb = int(15 * batch_size * n**2 * D * 32 / 1e6) + 1

    results = {}

    # 1. Naive autodiff (using raw function without custom_vjp)
    print("\n1. Naive autodiff...")
    naive_grad_fn = jax.jit(
        jax.grad(
            lambda p: jnp.real(
                jnp.sum(
                    _calculate_Gaussian_expectation_values_raw(p, pot_params, linear_params)
                )
            )
        )
    )
    t, grad_naive = time_gradient(naive_grad_fn, params)
    results["naive"] = t
    print(f"   Time: {t*1000:.3f} ms")

    # 2. Fast (fully vectorized custom vjp)
    print("\n2. Fast (vectorized custom VJP)...")
    fast_grad_fn = jax.jit(
        jax.grad(
            lambda p: jnp.real(
                jnp.sum(
                    calculate_Gaussian_expectation_values_fast(p, pot_params, linear_params)
                )
            )
        )
    )
    t, grad_fast = time_gradient(fast_grad_fn, params)
    results["fast"] = t
    print(f"   Time: {t*1000:.3f} ms")

    # 3. Batched (10 batches)
    print(f"\n3. Batched (10 batches, memory_limit={memory_limit_mb} MB)...")
    batched_grad_fn = jax.jit(
        jax.grad(
            lambda p: jnp.real(
                jnp.sum(
                    calculate_Gaussian_expectation_values_batched(
                        p, pot_params, linear_params, memory_limit_mb
                    )
                )
            )
        )
    )
    t, grad_batched = time_gradient(batched_grad_fn, params)
    results["batched"] = t
    print(f"   Time: {t*1000:.3f} ms")

    # Verify all gradients match
    print("Verification (max abs difference from naive):")
    print(f"   Fast:      {jnp.max(jnp.abs(grad_fast - grad_naive)):.2e}")
    print(f"   Batched:   {jnp.max(jnp.abs(grad_batched - grad_naive)):.2e}")

    return results


def main():
    print("Gradient Implementation Timing Comparison")
    print("=========================================")
    print("All times are for backward pass only (gradient computation)")

    # Problem sizes
    problems = [
        ("Tiny", 30, 55),
        ("Small", 50, 100),
        ("Medium", 60, 400),
    ]

    all_results = {}
    for name, n, m in problems:
        all_results[name] = run_benchmark(n, m)

    # Summary table
    print("" + "=" * 60)

    print("SUMMARY (times in ms)")
    print("=" * 70)
    print(f"{'Problem':<10} {'Naive':>10} {'Fast':>10} {'Batched':>10}")
    print("-" * 70)
    for name, n, m in problems:
        r = all_results[name]
        print(
            f"{name:<10} {r['naive']*1000:>10.3f} {r['fast']*1000:>10.3f} "
            f"{r['batched']*1000:>10.3f}"
        )

    # Speedup table
    print("" + "=" * 60)

    print("SPEEDUP vs Naive")
    print("=" * 70)
    print(f"{'Problem':<10} {'Fast':>10} {'Batched':>10}")
    print("-" * 70)
    for name, n, m in problems:
        r = all_results[name]
        print(
            f"{name:<10} {r['naive']/r['fast']:>10.2f}x " f"{r['naive']/r['batched']:>10.2f}x"
        )


if __name__ == "__main__":
    main()
