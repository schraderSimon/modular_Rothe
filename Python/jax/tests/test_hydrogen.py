"""Test hydrogen atom ground state energy and variance."""

import numpy as np
import pandas as pd
import pytest

import jax.numpy as jnp
from rothe.basis.linalg import calculate_variance, orthonormalize_matrices
from rothe.wavefunction import generalPotentialSolver


@pytest.fixture
def hydrogen_setup():
    """Set up hydrogen atom with Gaussian-expanded Coulomb potential."""
    csv_path = "gaussian_Coulomb/coeffs_mu=100_N=19.csv"
    df = pd.read_csv(csv_path)
    lincoeffs = df["linear"]
    width = df["nonlinear"]

    hydrogen_string = """
    dimension 3
    exponential
    """
    for lc, w in zip(lincoeffs, width):
        hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"

    n = 20
    random_widths = np.logspace(-2, 2, n, dtype=jnp.float64)
    n = len(random_widths)
    initial_params = jnp.zeros((n, 3, 4), dtype=jnp.float64)
    for i in range(n):
        initial_params = initial_params.at[i, :, 0].set(random_widths[i])

    hydrogen_wf = generalPotentialSolver(initial_params, initial_params, hydrogen_string)
    S = hydrogen_wf.calculate_S()
    S = S + 1e-8 * jnp.eye(S.shape[0], dtype=S.dtype)
    H = hydrogen_wf.calculate_H()
    H2 = hydrogen_wf.calculate_H2()

    result = orthonormalize_matrices(S, H, H2, return_Smin05=True)
    if len(result) == 5:
        eigvals, eigvecs_orth, H_orth, H2_orth, inverse = result
    else:
        eigvals, eigvecs_orth, H_orth, H2_orth = result
        inverse = None

    return eigvals, eigvecs_orth, H_orth, H2_orth


def test_hydrogen_ground_state_energy(hydrogen_setup):
    """Ground state energy should be close to -0.5 Hartree."""
    eigvals, _, _, _ = hydrogen_setup
    assert jnp.isclose(eigvals[0].real, -0.5, atol=1e-3), f"Expected E0 ~ -0.5, got {eigvals[0].real}"


def test_hydrogen_ground_state_variance(hydrogen_setup):
    """Ground state variance should be close to zero."""
    eigvals, eigvecs_orth, _, H2_orth = hydrogen_setup
    groundstate = eigvecs_orth[:, 0]
    variance = calculate_variance(groundstate, H2_orth, eigvals[0].real)
    assert variance >= 0
    assert variance < 1e-2, f"Variance too large: {variance}"
