"""
Centralized JAX/XLA configuration for the Rothe solver package.

This module is imported first by ``rothe/__init__.py`` so that every
downstream import sees the correct settings.  Do **not** set these
flags anywhere else in the codebase.
"""

import os

os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true --xla_gpu_autotune_level=4"

import jax  # noqa: E402 (must follow env-var setup)

jax.config.update("jax_enable_x64", True)
jax.config.update("jax_default_matmul_precision", "high")
