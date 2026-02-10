"""
Rothe solver for the time-dependent Schrödinger equation.

A JAX-based Gaussian-basis Rothe time-stepping method using implicit
variational steps optimized with BFGS.
"""

# Centralised JAX / XLA configuration — must be first import
from . import _config  # noqa: F401

# Public API
from .io import OutputConfig, load_rothe_file, save_rothe_state
from .solver import RotheSolver, setUpRotheErrorAndGradient_jit
from .wavefunction import (
    ND_potentials,
    generalPotentialSolver,
    propagate_kinetic_analytical,
)

__all__ = [
    "RotheSolver",
    "setUpRotheErrorAndGradient_jit",
    "OutputConfig",
    "save_rothe_state",
    "load_rothe_file",
    "ND_potentials",
    "generalPotentialSolver",
    "propagate_kinetic_analytical",
]
