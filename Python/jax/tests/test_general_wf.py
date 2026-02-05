import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

"""
Tests for general_wf.py functions
"""

import numpy as np

import jax.numpy as jnp
from libraries.general_wf import (
    calculate_Gaussian_expectation_values,
    calculate_squared_Gaussian_potential,
    dtype_complex,
    dtype_real,
    generalPotentialSolver,
)
from libraries.quadrature_utils import (
    eval_gaussian_basis,
    eval_gaussian_potential,
    eval_polynomial,
    quadrature_cross_terms,
    quadrature_gaussian_kinetic_cross,
)


def test_squared_gaussian_potential():
    """
    Test calculate_squared_Gaussian_potential by evaluating both the original
    Gaussians on a grid and squaring them, then comparing with the function's output.
    """
    np.random.seed(42)

    # Create 4 random Gaussians in 2D
    m, D = 4, 2
    exponential_params = jnp.zeros((m, D, 4), dtype=dtype_real)
    linear_params = jnp.zeros(m, dtype=dtype_complex)

    # Random parameters
    for k in range(m):
        a_k = np.random.uniform(0.5, 2.0, D)  # width parameters
        mu_k = np.random.uniform(-2.0, 2.0, D)  # centers
        p_k = np.zeros(D)  # no momentum
        c_k = np.random.uniform(0.5, 1.5) + 1j * np.random.uniform(-0.5, 0.5)

        exponential_params = exponential_params.at[k, :, 0].set(a_k)
        exponential_params = exponential_params.at[k, :, 1].set(0)  # b=0
        exponential_params = exponential_params.at[k, :, 2].set(mu_k)
        exponential_params = exponential_params.at[k, :, 3].set(p_k)
        linear_params = linear_params.at[k].set(c_k)

    # Compute squared potential using the function
    sq_params, sq_linear = calculate_squared_Gaussian_potential(
        exponential_params, linear_params
    )

    # Create evaluation grid
    x_range = np.linspace(-5, 5, 30)
    y_range = np.linspace(-5, 5, 30)
    xx, yy = np.meshgrid(x_range, y_range)

    # Evaluate original Gaussians on grid
    def eval_gaussian(x, y, a, b, mu, p, c):
        """Evaluate single Gaussian at (x, y)"""
        # a, b, mu, p are arrays of shape (D,) - need to sum over dimensions
        alpha = a**2 + 1j * b  # shape (D,)
        # Sum exponents over dimensions
        exponent = -alpha[0] * (x - mu[0]) ** 2 - alpha[1] * (y - mu[1]) ** 2
        return c * np.exp(exponent)

    # Sum over all original Gaussians squared
    V_squared_direct = np.zeros_like(xx, dtype=dtype_complex)
    for i in range(m):
        for j in range(m):
            a_i = exponential_params[i, :, 0]
            b_i = exponential_params[i, :, 1]
            mu_i = exponential_params[i, :, 2]
            p_i = exponential_params[i, :, 3]
            c_i = linear_params[i]

            a_j = exponential_params[j, :, 0]
            b_j = exponential_params[j, :, 1]
            mu_j = exponential_params[j, :, 2]
            p_j = exponential_params[j, :, 3]
            c_j = linear_params[j]

            gi = eval_gaussian(xx, yy, a_i, b_i, mu_i, p_i, c_i)
            gj = eval_gaussian(xx, yy, a_j, b_j, mu_j, p_j, c_j)
            V_squared_direct += gi * gj

    # Evaluate using squared Gaussian params
    V_squared_reconstructed = np.zeros_like(xx, dtype=dtype_complex)
    n_sq = sq_params.shape[0]
    for k in range(n_sq):
        a_k = sq_params[k, :, 0]
        b_k = sq_params[k, :, 1]
        mu_k = sq_params[k, :, 2]
        p_k = sq_params[k, :, 3]
        c_k = sq_linear[k]

        gk = eval_gaussian(xx, yy, a_k, b_k, mu_k, p_k, c_k)
        V_squared_reconstructed += gk

    # Compare
    max_error = np.max(np.abs(V_squared_direct - V_squared_reconstructed))
    rel_error = max_error / (np.max(np.abs(V_squared_direct)) + 1e-15)

    print("\n" + "=" * 60)
    print("TEST: calculate_squared_Gaussian_potential")
    print("=" * 60)
    print(f"Number of original Gaussians: {m}")
    print(f"Dimension: {D}")
    print(f"Number of squared Gaussians: {n_sq} (expected {m*(m+1)//2})")
    print(f"Max absolute error: {max_error:.2e}")
    print(f"Relative error: {rel_error:.2e}")

    if rel_error < 1e-10:
        print("✓ PASSED: Squared Gaussians reconstruct correctly!")
        return True
    else:
        print("✗ FAILED: Large discrepancy in reconstruction")
        print(f"Sample direct values:\\n{V_squared_direct[10:15, 10:15]}")
        print(f"Sample reconstructed values:\\n{V_squared_reconstructed[10:15, 10:15]}")
        return False


def test_gaussian_expectation_values():
    """Test basic Gaussian expectation value calculation"""
    n = 2
    D = 4
    params_init = jnp.zeros((n, D, 4), dtype=dtype_real)
    params_init = params_init.at[:, :, 0].set(1 / jnp.sqrt(2))  # Width parameters
    params_init = params_init.at[0, :, 2].set(2.0)  # Set mu to (2,...,2)

    example_string = """
    dimension 2
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x0x0x1: 0.111803
    x1x1x1: -0.03726766666
    x0x0x0x0: 0.0007812444255625
    x1x1x1x1: 0.0007812444255625
    x0x0x1x1: 0.001562488851125
    """
    for i in range(1, n):
        params_init = params_init.at[i, :, 2].set(
            2 + np.random.uniform(-0.5, 0.5, (D,))
        )  # Move to the right

    osc = generalPotentialSolver(params_init, example_string)
    S = osc.calculate_S()
    gaussian = jnp.zeros((2, D, 4), dtype=dtype_real)
    gaussian = gaussian.at[0, :, 0].set(1e-11)  # a
    gaussian = gaussian.at[1, :, 0].set(1e-11)  # b
    linear_params = jnp.ones(gaussian.shape[0], dtype=dtype_complex)
    gev = calculate_Gaussian_expectation_values(params_init, gaussian, linear_params)

    print("\\n" + "=" * 60)
    print("TEST: Gaussian expectation values")
    print("=" * 60)
    print("S=", S)
    print("Gaussian expectation values:", gev)
    print("✓ Test completed (visual inspection required)")


def test_moments():
    """Test moment calculation"""
    n = 3
    D = 2
    params_init = jnp.zeros((n, D, 4), dtype=dtype_real)
    params_init = params_init.at[0, 0, :].set([1 / jnp.sqrt(2), 0, 1, 0])
    params_init = params_init.at[0, 1, :].set([1 / jnp.sqrt(2), 0, 0, 0])
    params_init = params_init.at[1, 0, :].set([1 / jnp.sqrt(2), 0, 1, 0])  # x expectation is 1
    params_init = params_init.at[1, 1, :].set([1 / jnp.sqrt(2), 0, -10, 0])
    params_init = params_init.at[2, 0, :].set([1 / jnp.sqrt(2), 0, 0, 0])  # x expectation is 1
    params_init = params_init.at[2, 1, :].set([1 / jnp.sqrt(2), 0, 0, 0])

    example_string = """
    dimension 2
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    """

    osc = generalPotentialSolver(params_init, example_string)
    S = osc.calculate_S()
    expectations = osc.moments

    print("" + "=" * 60)
    print("TEST: Moment calculation")
    print("=" * 60)
    print("Parameters:")
    print(params_init)
    print("First moment (x expectation):")
    print(expectations[:, :, 1, 1])
    print("✓ Test completed (visual inspection required)")


def test_polynomial_gaussian_cross_terms_quadrature():
    """Numerically validate polynomial–Gaussian cross terms in 2D via quadrature."""

    # Complex basis (nonzero b,p) with real potential (b=p=0)
    n, D = 2, 2
    params = jnp.zeros((n, D, 4), dtype=dtype_real)
    params = params.at[:, :, 0].set(0.8)  # width a
    params = params.at[:, :, 1].set(jnp.array([[0.3, -0.25], [-0.15, 0.2]]))
    params = params.at[0, :, 2].set(jnp.array([0.0, -0.5]))
    params = params.at[1, :, 2].set(jnp.array([0.7, 0.6]))
    params = params.at[0, :, 3].set(jnp.array([0.4, -0.2]))
    params = params.at[1, :, 3].set(jnp.array([-0.3, 0.35]))

    polynomial_string = """
    dimension 2
    polynomial
    x0: 0.3
    x1x1: -0.2
    x0x1: 0.12
    x0x0x1: -0.05
    exponential
    0.9, 0.6, [ -0.4, 0.8 ]
    0.4, 0.7, [ 1.1, -0.3 ]
    0.3, 0.9, [ -1.2, -0.6 ]
    """

    # Parse potentials directly from string (no manual assignment)
    solver = generalPotentialSolver(params, polynomial_string)
    exp_params = solver.exponential_params
    lin_params = solver.linear_exponential_params

    S = solver.calculate_S()
    analytic_cross = solver.calculate_polynomial_gaussian_cross_terms(S)

    # Quadrature grid
    L = 6.0
    N = 181  # fine grid for reasonable accuracy
    xs = np.linspace(-L, L, N)
    ys = np.linspace(-L, L, N)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    X, Y = np.meshgrid(xs, ys, indexing="ij")

    numeric_cross = quadrature_cross_terms(
        params=params,
        poly_dict=solver.polynomial,
        exp_params=exp_params,
        lin_params=lin_params,
        xs=xs,
        ys=ys,
    )
    print("Numeric cross terms:")
    print(numeric_cross)
    print("Analytic cross terms:")
    print(np.array(analytic_cross))
    max_err = np.max(np.abs(numeric_cross - np.array(analytic_cross)))
    rel_err = max_err / (np.max(np.abs(numeric_cross)) + 1e-14)

    print("" + "=" * 60)

    print("TEST: polynomial–Gaussian cross terms (quadrature vs analytic)")
    print("=" * 60)
    print(f"Max abs error: {max_err:.3e}")
    print(f"Rel error:     {rel_err:.3e}")

    assert rel_err < 5e-3, "Cross term deviates from quadrature reference"


def test_gaussian_kinetic_cross_terms_quadrature():
    """Quadrature check for Gaussian–kinetic cross term {T,V_gauss} in 2D."""

    n, D = 2, 2
    params = jnp.zeros((n, D, 4), dtype=dtype_real)
    params = params.at[:, :, 0].set(0.9)
    params = params.at[:, :, 1].set(jnp.array([[0.2, -0.15], [-0.1, 0.18]]))
    params = params.at[0, :, 2].set(jnp.array([-0.2, 0.4]))
    params = params.at[1, :, 2].set(jnp.array([0.9, -0.3]))
    params = params.at[0, :, 3].set(jnp.array([0.35, -0.25]))
    params = params.at[1, :, 3].set(jnp.array([-0.28, 0.42]))

    polynomial_string = """
    dimension 2
    polynomial
    x0: 0.1
    x1: -0.05
    exponential
    1.1, 0.7, [ -0.5, 0.9 ]
    0.6, 0.8, [ 1.0, -0.4 ]
    """

    solver = generalPotentialSolver(params, polynomial_string)
    exp_params = solver.exponential_params
    lin_params = solver.linear_exponential_params

    S = solver.calculate_S()
    analytic = solver.calculate_gaussian_kinetic_cross_terms(S)

    L = 6.0
    N = 181
    xs = np.linspace(-L, L, N)
    ys = np.linspace(-L, L, N)

    numeric = quadrature_gaussian_kinetic_cross(params, exp_params, lin_params, xs, ys)

    max_err = np.max(np.abs(numeric - np.array(analytic)))
    rel_err = max_err / (np.max(np.abs(numeric)) + 1e-14)

    print("" + "=" * 60)

    print("TEST: Gaussian–kinetic cross term (quadrature vs analytic)")
    print("=" * 60)
    print(f"Max abs error: {max_err:.3e}")
    print(f"Rel error:     {rel_err:.3e}")

    assert rel_err < 7e-3, "Gaussian–kinetic cross term deviates from quadrature reference"


if __name__ == "__main__":
    # Run all tests
    test_squared_gaussian_potential()
    test_gaussian_expectation_values()
    test_moments()
    test_polynomial_gaussian_cross_terms_quadrature()
    test_gaussian_kinetic_cross_terms_quadrature()
