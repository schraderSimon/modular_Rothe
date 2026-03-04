"""
Henon-Heiles system helpers.

Provides potential strings for 1D, 2D, 4D, and 6D Henon-Heiles chains,
plus parameter and coefficient initialization routines.
"""

import numpy as np

import jax.numpy as jnp

# ---------------------------------------------------------------------------
# Pre-built potential strings for standard Henon-Heiles chain dimensions
# ---------------------------------------------------------------------------


def make_henon_heiles_string(dim: int) -> str:
    fourth_order_term = 0.001562488851125
    val = fourth_order_term / dim  #
    string = f"""dimension {dim}\npolynomial\n"""
    for i in range(dim):
        string += f"x{i}x{i}: 0.5\n"  # Quadratic terms
    for i in range(dim - 1):
        string += f"x{i}x{i}x{i + 1}: 0.111803\n"  # Cubic coupling terms
    for i in range(1, dim):
        string += f"x{i}x{i}x{i}: -0.03726766666\n"  # Cubic self-interaction terms
    for i in range(dim):
        string += f"x{i}x{i}x{i}x{i}: {val}\n"  # Quartic safety terms to avoid unboundedness
    print(string)
    return string


def initialize_henon_heiles_params(n: int, D: int, seed: int = 42):
    """
    Initialize parameters for Henon-Heiles system.
    Gaussians will be placed symmetrically around (2, 2, ..., 2) with some random perturbation,
    with one Gaussian exactly at (2, 2, ..., 2).
    Gaussians will have width 1/√2 in all dimensions and zero momentum.
    Parameters:
        n: number of Gaussians
        D: dimension
        seed: random seed

    Returns:
        params_init: array of shape (n, D, 4)
    """
    np.random.seed(seed)
    params_init = jnp.zeros((n, D, 4), dtype=jnp.float64)
    params_init = params_init.at[:, :, 0].set(1 / jnp.sqrt(2))
    params_init = params_init.at[0, :, 2].set(2.0)
    for i in range(1, n):
        params_init = params_init.at[i, :, 2].set(2 + np.random.uniform(-0.5, 0.5, (D,)))
    return params_init
