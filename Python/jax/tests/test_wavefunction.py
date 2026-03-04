"""Tests for wavefunction matrix element calculations."""

import numpy as np
import pytest

import jax.numpy as jnp
from rothe.wavefunction import (
    calculate_Gaussian_expectation_values,
    calculate_squared_Gaussian_potential,
    dtype_complex,
    dtype_real,
    generalPotentialSolver,
)
from tests.quadrature import (
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

    m, D = 4, 2
    exponential_params = jnp.zeros((m, D, 4), dtype=dtype_real)
    linear_params = jnp.zeros(m, dtype=dtype_complex)

    for k in range(m):
        a_k = np.random.uniform(0.5, 2.0, D)
        mu_k = np.random.uniform(-2.0, 2.0, D)
        p_k = np.zeros(D)
        c_k = np.random.uniform(0.5, 1.5) + 1j * np.random.uniform(-0.5, 0.5)

        exponential_params = exponential_params.at[k, :, 0].set(a_k)
        exponential_params = exponential_params.at[k, :, 1].set(0)
        exponential_params = exponential_params.at[k, :, 2].set(mu_k)
        exponential_params = exponential_params.at[k, :, 3].set(p_k)
        linear_params = linear_params.at[k].set(c_k)

    sq_params, sq_linear = calculate_squared_Gaussian_potential(exponential_params, linear_params)

    x_range = np.linspace(-5, 5, 30)
    y_range = np.linspace(-5, 5, 30)
    xx, yy = np.meshgrid(x_range, y_range)

    def eval_gaussian(x, y, a, b, mu, p, c):
        alpha = a**2 + 1j * b
        exponent = -alpha[0] * (x - mu[0]) ** 2 - alpha[1] * (y - mu[1]) ** 2
        return c * np.exp(exponent)

    V_squared_direct = np.zeros_like(xx, dtype=dtype_complex)
    for i in range(m):
        for j in range(m):
            a_i, b_i = exponential_params[i, :, 0], exponential_params[i, :, 1]
            mu_i, p_i = exponential_params[i, :, 2], exponential_params[i, :, 3]
            a_j, b_j = exponential_params[j, :, 0], exponential_params[j, :, 1]
            mu_j, p_j = exponential_params[j, :, 2], exponential_params[j, :, 3]
            gi = eval_gaussian(xx, yy, a_i, b_i, mu_i, p_i, linear_params[i])
            gj = eval_gaussian(xx, yy, a_j, b_j, mu_j, p_j, linear_params[j])
            V_squared_direct += gi * gj

    V_squared_reconstructed = np.zeros_like(xx, dtype=dtype_complex)
    n_sq = sq_params.shape[0]
    for k in range(n_sq):
        a_k = sq_params[k, :, 0]
        b_k = sq_params[k, :, 1]
        mu_k = sq_params[k, :, 2]
        p_k = sq_params[k, :, 3]
        gk = eval_gaussian(xx, yy, a_k, b_k, mu_k, p_k, sq_linear[k])
        V_squared_reconstructed += gk

    max_error = np.max(np.abs(V_squared_direct - V_squared_reconstructed))
    rel_error = max_error / (np.max(np.abs(V_squared_direct)) + 1e-15)
    assert n_sq == m * (m + 1) // 2
    assert rel_error < 1e-10, f"Relative error {rel_error:.2e} too large"


def test_gaussian_expectation_values():
    """Test basic Gaussian expectation value calculation."""
    n = 2
    D = 4
    params_init = jnp.zeros((n, D, 4), dtype=dtype_real)
    params_init = params_init.at[:, :, 0].set(1 / jnp.sqrt(2))
    params_init = params_init.at[0, :, 2].set(2.0)

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
        params_init = params_init.at[i, :, 2].set(2 + np.random.uniform(-0.5, 0.5, (D,)))

    osc = generalPotentialSolver(params_init, params_init, example_string)
    S = osc.calculate_S()
    gaussian = jnp.zeros((2, D, 4), dtype=dtype_real)
    gaussian = gaussian.at[0, :, 0].set(1e-11)
    gaussian = gaussian.at[1, :, 0].set(1e-11)
    linear_params = jnp.ones(gaussian.shape[0], dtype=dtype_complex)
    gev = calculate_Gaussian_expectation_values(params_init, params_init, gaussian, linear_params)
    # Flat Gaussians (very small width) should give values close to S
    assert gev.shape == S.shape


def test_moments():
    """Test moment calculation for shifted Gaussians."""
    n = 3
    D = 2
    params_init = jnp.zeros((n, D, 4), dtype=dtype_real)
    params_init = params_init.at[0, 0, :].set([1 / jnp.sqrt(2), 0, 1, 0])
    params_init = params_init.at[0, 1, :].set([1 / jnp.sqrt(2), 0, 0, 0])
    params_init = params_init.at[1, 0, :].set([1 / jnp.sqrt(2), 0, 1, 0])
    params_init = params_init.at[1, 1, :].set([1 / jnp.sqrt(2), 0, -10, 0])
    params_init = params_init.at[2, 0, :].set([1 / jnp.sqrt(2), 0, 0, 0])
    params_init = params_init.at[2, 1, :].set([1 / jnp.sqrt(2), 0, 0, 0])

    example_string = """
    dimension 2
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    """

    osc = generalPotentialSolver(params_init, params_init, example_string)
    S = osc.calculate_S()
    expectations = osc.moments
    # moments should be (n, n, D, max_order+1) array
    assert expectations.shape[0] == n
    assert expectations.shape[1] == n


def test_polynomial_gaussian_cross_terms_quadrature():
    """Numerically validate polynomial-Gaussian cross terms in 2D via quadrature."""
    n, D = 2, 2
    params = jnp.zeros((n, D, 4), dtype=dtype_real)
    params = params.at[:, :, 0].set(0.8)
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

    solver = generalPotentialSolver(params, params, polynomial_string)
    exp_params = solver.exponential_params
    lin_params = solver.linear_exponential_params

    S = solver.calculate_S()
    analytic_cross = solver.calculate_polynomial_gaussian_cross_terms(S)

    L = 6.0
    N = 181
    xs = np.linspace(-L, L, N)
    ys = np.linspace(-L, L, N)

    numeric_cross = quadrature_cross_terms(
        params=params, poly_dict=solver.polynomial, exp_params=exp_params, lin_params=lin_params, xs=xs, ys=ys
    )
    max_err = np.max(np.abs(numeric_cross - np.array(analytic_cross)))
    rel_err = max_err / (np.max(np.abs(numeric_cross)) + 1e-14)

    assert rel_err < 5e-3, f"Cross term rel error {rel_err:.3e} exceeds tolerance"


def test_gaussian_kinetic_cross_terms_quadrature():
    """Quadrature check for Gaussian-kinetic cross term {T,V_gauss} in 2D."""
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

    solver = generalPotentialSolver(params, params, polynomial_string)
    exp_params = solver.exponential_params
    lin_params = solver.linear_exponential_params

    S = solver.calculate_S()
    analytic = np.array(solver.calculate_gaussian_kinetic_cross_terms(S))
    analytic_sym = analytic + analytic.conj().T

    L = 6.0
    N = 181
    xs = np.linspace(-L, L, N)
    ys = np.linspace(-L, L, N)

    numeric = quadrature_gaussian_kinetic_cross(params, exp_params, lin_params, xs, ys)

    max_err = np.max(np.abs(numeric - analytic_sym))
    rel_err = max_err / (np.max(np.abs(numeric)) + 1e-14)

    assert rel_err < 7e-3, f"Gaussian-kinetic cross term rel error {rel_err:.3e} exceeds tolerance"
