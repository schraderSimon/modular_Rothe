from math import comb

import numpy as np

import jax.numpy as jnp
from general_wf import generalPotentialSolver
from utils import _symh, calculate_variance, eigh_canonical, orthonormalize_matrices


def create_harmonic_oscillator_string(D, omega=1.0):
    lines = [f"    dimension {D}", "    polynomial"]
    coeff = 0.5 * omega**2
    for i in range(D):
        lines.append(f"    x{i}x{i}: {coeff}")
    return "\n".join(lines) + "\n"


def groundstate_energy_nD(params, polynomial_string, eps=1e-12):
    """
    Calculate ground state energy for n-dimensional system.
    """
    solver = generalPotentialSolver(params, polynomial_string)
    S, H, H2 = solver.calculate_SHH2()

    # Canonical orthogonalization and diagonalization
    eigvals, eigvecs = eigh_canonical(S, H, eps=eps)

    return jnp.real(eigvals[0])


def degeneracy(n, D):

    return comb(n + D - 1, D - 1)


def initialize_gaussian_parameters(num_gaussians, D, seed=42):
    """
    Create a set of randomly perturbed Gaussian parameters.
    """
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


def print_energy_level_summary(level, expected, deg, D):
    print(f"      Level {level} (expected E={expected:.1f}, degeneracy={deg}):")


def test_groundstate_nD(Dmax):
    """Test that the ground state has energy 0.5*D for various dimensions."""
    print("Testing ground state energy for exact ground state...")
    Dvals = range(1, Dmax + 1)
    for D in Dvals:
        # Create exact ground state: single Gaussian centered at origin
        params = jnp.zeros((1, D, 4))
        params = params.at[:, :, 0].set(1 / jnp.sqrt(2))  # a = 1/sqrt(2)

        polynomial_string = create_harmonic_oscillator_string(D, omega=1.0)
        solver = generalPotentialSolver(params, polynomial_string)
        S, H, H2 = solver.calculate_SHH2()

        energy = H[0, 0]
        expected_energy = 0.5 * D

        print(f"  D={D}: Energy = {energy:.6f}, Expected = {expected_energy:.6f}")
        assert jnp.isclose(
            energy, expected_energy, rtol=1e-10
        ), f"Ground state energy mismatch for D={D}: got {energy}, expected {expected_energy}"

        # Also check that variance is zero for exact ground state
        variance = H2[0, 0] - energy**2
        print(f"        Variance = {variance:.2e}")
        assert jnp.isclose(
            variance, 0.0, atol=1e-10
        ), f"Variance should be zero for exact ground state, got {variance}"

        # Check using the groundstate_energy function
        gs_energy = groundstate_energy_nD(params, polynomial_string)
        assert jnp.isclose(gs_energy, expected_energy, rtol=1e-10)

    print("✓ Ground state energy test passed!\n")


def _test_dimension(D, num_gaussians, num_levels):
    """
    Test a specific dimension with the given basis size,
    based on a random initialization of Gaussian parameters.
    """
    print(f"\n  Testing D={D} with {num_gaussians} Gaussians:")

    polynomial_string = create_harmonic_oscillator_string(D, omega=1.0)

    params = initialize_gaussian_parameters(num_gaussians, D, seed=42)

    solver = generalPotentialSolver(params, polynomial_string)
    S, H, H2 = solver.calculate_SHH2()

    eigvals, eigvecs_orth, H_orth, H2_orth = orthonormalize_matrices(S, H, H2)

    print(f"    First {num_levels} energy levels:")

    expected_energies = jnp.arange(num_levels) + 0.5 * D

    # Check each energy level
    idx = 0
    for level in range(num_levels):
        expected = expected_energies[level]
        deg = degeneracy(level, D)

        print_energy_level_summary(level, expected, deg, D)

        # Check all degenerate states
        for i in range(deg):
            if idx + i >= len(eigvals):
                print(f"        WARNING: Not enough eigenvalues computed!")
                break

            computed = eigvals[idx + i]
            error = abs(computed - expected)

            # Calculate variance
            state_vec = eigvecs_orth[:, idx + i]
            var = calculate_variance(state_vec, H2_orth, computed)

            print(f"        State {i}: E={computed:.6f}, error={error:.2e}, var={var:.2e}")

            # Tolerance depends on basis size and level
            rtol = 1e-2
            assert jnp.isclose(computed, expected, rtol=rtol), (
                f"Energy mismatch for D={D}, level={level}, state={i}: "
                f"got {computed}, expected {expected}"
            )

        idx += deg

    print(f"  ✓ D={D} test passed!")


def test_perturbed_gaussians_nD(Dmax):
    """
    Test that using many perturbed Gaussians recovers correct low-lying energies.
    For a D-dimensional harmonic oscillator, the energy levels are:
    E_n = (n + D/2) for n = 0, 1, 2, ...
    where n is the total quantum number (sum of individual mode quantum numbers).

    Degeneracy of level n is C(n+D-1, D-1) = (n+D-1)! / (n! * (D-1)!)
    """
    print("Testing energy levels with perturbed Gaussians...")

    # Test parameters for different dimensions
    for D in range(1, Dmax + 1):
        num_gaussians = 20 * D * (D + 1)  # Increase basis size with dimension
        num_levels = 3
        _test_dimension(D, num_gaussians, num_levels)

    print("\n✓ All perturbed Gaussian tests passed!")


if __name__ == "__main__":
    test_groundstate_nD(5)
    test_perturbed_gaussians_nD(5)
    print("\n" + "=" * 60)
    print("All n-dimensional harmonic oscillator tests passed!")
    print("=" * 60)
