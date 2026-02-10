"""
Gaussian basis integral machinery.

Subpackage containing:
- ``gaussian_potential``: Batched <g_i|V|g_j> for Gaussian potentials (custom VJP).
- ``hessian``: SPD Hessian approximation from overlap moments.
- ``linalg``: Canonical orthogonalization and related utilities.
"""

from .gaussian_potential import (
    calculate_Gaussian_expectation_values_batched,
    calculate_squared_Gaussian_potential,
    get_exponent_subtraction_terms,
    get_gaussian_pair_terms,
    parse_exponential_params,
)
from .linalg import calculate_variance, eigh_canonical, orthonormalize_matrices, symh

__all__ = [
    "calculate_Gaussian_expectation_values_batched",
    "calculate_squared_Gaussian_potential",
    "get_exponent_subtraction_terms",
    "get_gaussian_pair_terms",
    "parse_exponential_params",
    "eigh_canonical",
    "orthonormalize_matrices",
    "calculate_variance",
    "symh",
]
