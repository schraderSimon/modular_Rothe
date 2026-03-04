"""Tests for n-dimensional harmonic oscillator energy levels."""

from math import comb

import numpy as np
import pytest

import jax.numpy as jnp
from rothe.basis.linalg import calculate_variance, eigh_canonical, orthonormalize_matrices
from rothe.wavefunction import generalPotentialSolver


def create_harmonic_oscillator_string(D, omega=1.0):
    lines = [f"    dimension {D}", "    polynomial"]
    coeff = 0.5 * omega**2
    for i in range(D):
        lines.append(f"    x{i}x{i}: {coeff}")
    return "\n".join(lines) + "\n"


def groundstate_energy_nD(params, polynomial_string, eps=1e-12):
    """Calculate ground state energy for n-dimensional system."""
    solver = generalPotentialSolver(params, params, polynomial_string)
    S, H, H2 = solver.calculate_SHH2()
    eigvals, eigvecs = eigh_canonical(S, H, eps=eps)
    return jnp.real(eigvals[0])


def degeneracy(n, D):
    return comb(n + D - 1, D - 1)


def initialize_gaussian_parameters(num_gaussians, D, seed=42):
    """Create a set of randomly perturbed Gaussian parameters."""
    np.random.seed(seed)
    params = jnp.zeros((num_gaussians, D, 4))

    base_a = 1 / jnp.sqrt(2)
    a_values = base_a * (1 + np.random.uniform(-0.15, 0.15, (num_gaussians, D)))
    params = params.at[:, :, 0].set(a_values)

    b_values = np.random.uniform(-0.1, 0.1, (num_gaussians, D))
    params = params.at[:, :, 1].set(b_values)

    mu_values = np.random.uniform(-1.0, 1.0, (num_gaussians, D))
    params = params.at[:, :, 2].set(mu_values)

    p_values = np.random.uniform(-0.5, 0.5, (num_gaussians, D))
    params = params.at[:, :, 3].set(p_values)

    return params


@pytest.mark.parametrize("D", range(1, 6))
def test_groundstate_nD(D):
    """Test that the ground state has energy 0.5*D."""
    params = jnp.zeros((1, D, 4))
    params = params.at[:, :, 0].set(1 / jnp.sqrt(2))

    polynomial_string = create_harmonic_oscillator_string(D, omega=1.0)
    solver = generalPotentialSolver(params, params, polynomial_string)
    S, H, H2 = solver.calculate_SHH2()

    energy = H[0, 0]
    expected_energy = 0.5 * D

    assert jnp.isclose(
        energy, expected_energy, rtol=1e-10
    ), f"Ground state energy mismatch for D={D}: got {energy}, expected {expected_energy}"

    variance = H2[0, 0] - energy**2
    assert jnp.isclose(variance, 0.0, atol=1e-10), f"Variance should be zero for exact ground state, got {variance}"

    gs_energy = groundstate_energy_nD(params, polynomial_string)
    assert jnp.isclose(gs_energy, expected_energy, rtol=1e-10)


@pytest.mark.parametrize("D", range(1, 6))
def test_perturbed_gaussians_nD(D):
    """Test that many perturbed Gaussians recover correct low-lying energies."""
    num_gaussians = 20 * D * (D + 1)
    num_levels = 3

    polynomial_string = create_harmonic_oscillator_string(D, omega=1.0)
    params = initialize_gaussian_parameters(num_gaussians, D, seed=42)

    solver = generalPotentialSolver(params, params, polynomial_string)
    S, H, H2 = solver.calculate_SHH2()

    eigvals, eigvecs_orth, H_orth, H2_orth = orthonormalize_matrices(S, H, H2)

    expected_energies = jnp.arange(num_levels) + 0.5 * D

    idx = 0
    for level in range(num_levels):
        expected = expected_energies[level]
        deg = degeneracy(level, D)

        for i in range(deg):
            if idx + i >= len(eigvals):
                break

            computed = eigvals[idx + i]
            state_vec = eigvecs_orth[:, idx + i]
            var = calculate_variance(state_vec, H2_orth, computed)

            assert jnp.isclose(computed, expected, rtol=1e-2), (
                f"Energy mismatch for D={D}, level={level}, state={i}: " f"got {computed}, expected {expected}"
            )

        idx += deg
