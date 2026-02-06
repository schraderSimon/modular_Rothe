import argparse
import os

os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true --xla_gpu_autotune_level=4"

import henon_heiles_strings
from jax import numpy as jnp
from libraries.file_handling import OutputConfig
from libraries.jax_Rothe import RotheSolver
from libraries.read_string import read_string
from libraries.system_helpers import (
    initialize_henon_heiles_coeffs,
    initialize_henon_heiles_params,
    resume_or_initialize_rothe,
    set_up_functions,
)


def parse_args():
    """Parse command-line arguments for Henon-Heiles simulation."""
    parser = argparse.ArgumentParser(
        description="Simulate Henon-Heiles system using the Rothe method."
    )
    parser.add_argument(
        "--dim",
        type=int,
        default=2,
        choices=[1, 2, 4, 6],
        help="Dimension of the system (1, 2, 4, or 6)",
    )
    parser.add_argument(
        "--splitting_type",
        type=str,
        default="none",
        help="Splitting type for Rothe propagation",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=0.01,
        help="Time step",
    )
    parser.add_argument(
        "--num_gaussians",
        type=int,
        default=1,
        help="Number of Gaussians",
    )
    return parser.parse_args()


# Load/set up parameters
args = parse_args()
splitting_type = args.splitting_type
dim = args.dim
dt = args.dt
n = args.num_gaussians
nsteps = int(100 / dt)

# System configuration
example_string = henon_heiles_strings.henonheiles_strings[f"{dim}D"]
polynomial, exponential = read_string(example_string)

# Extract dimension and number of Gaussians from polynomial
key, value = next(iter(polynomial.items()))
D = len(key)

# Initialize parameters
params_init = initialize_henon_heiles_params(n=n, D=D, seed=42)

# Set up functions and solver
rothe_error, rothe_vg_jit, SHH2 = set_up_functions(
    params_init, example_string, splitting_type=splitting_type
)

# Set up wavefunction coefficients
coeffs_init = initialize_henon_heiles_coeffs(
    n=n, params_init=params_init, potential_string=example_string
)

# Calculate initial error
initial_error = rothe_error(params_init, params_init, coeffs_init, SHH2, 0, dt)

# Set up Rothe solver
rothesolver = RotheSolver(
    SHH2,
    dt,
    0,
    params_init,
    coeffs_init,
    rothe_grad_jax=rothe_vg_jit,
    rothe_nograd=rothe_error,
    splitting_type=splitting_type,
    output_config=OutputConfig(
        name="HenonHeiles_%dD_%d" % (D, n),
        polynomial_string=example_string,
        out_dir="./wave_function_data",
    ),
)

# Propagation setup
nsteps = resume_or_initialize_rothe(
    rothesolver,
    nsteps,
    rothesolver.name,
    splitting_type,
    example_string,
    dt,
    initial_error=initial_error,
    params=params_init,
    coeffs=coeffs_init,
)

rothesolver.propagate(nsteps, maxiter=1000)
