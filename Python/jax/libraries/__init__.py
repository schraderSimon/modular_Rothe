"""
Core libraries for the modular Rothe solver.

This package contains the main computational modules for the JAX-based
Rothe time-stepping method for solving the time-dependent Schrödinger equation.
"""

# Just make the modules available, don't import specific functions
# This avoids issues with missing functions or circular imports

__all__ = [
    "general_wf",
    "jax_Rothe",
    "gaussian_potential_helpers",
    "file_handling",
    "optimization_helpers",
    "helpers",
    "polynomial_kinetic_vjp",
    "quadrature_utils",
    "read_string",
    "utils",
    "calculate_Hessian_coefficients",
]
