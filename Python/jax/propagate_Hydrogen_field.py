import argparse

import numpy as np

from initial_state_calculations.groundstate_from_imaginary_data import (
    read_best_coefficients,
)
from jax import numpy as jnp
from libraries.file_handling import OutputConfig
from libraries.general_wf import generalPotentialSolver
from libraries.jax_Rothe import RotheSolver
from libraries.system_helpers import (
    ExternalFieldParams,
    make_hydrogen_string,
    resume_or_initialize_rothe,
    set_up_functions,
    set_up_hydrogen_string_name,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Propagate the hydrogen atom wave function using the Rothe method."
    )
    parser.add_argument("--num_gauss_potential", type=int, default=21)
    parser.add_argument("--num_gauss_wavefunction", type=int)
    parser.add_argument("--dt", type=complex, help="time step", default=0.2)
    parser.add_argument("--rothe_epsilon", type=float, default=0.5)
    parser.add_argument("--maxiter", type=int, default=50)
    parser.add_argument(
        "--t_start",
        type=float,
        default=None,
        help="Time to resume propagation from (resumes from checkpoint at or before this time)",
    )
    return parser.parse_args()


# Load/set up parameters
args = parse_args()
ngp = num_gauss_potential = args.num_gauss_potential
ngwf = num_gauss_wavefunction = args.num_gauss_wavefunction
rothe_epsilon = args.rothe_epsilon  # Unused right now
maxiter = args.maxiter  # Unused right now
dt = args.dt
t_start = args.t_start
external_field_params = ExternalFieldParams()  # Standard parameters
string_name = set_up_hydrogen_string_name(external_field_params, ngp, ngwf, dt)
tfinal = external_field_params.t_final
nsteps = int(jnp.ceil(tfinal / jnp.abs(dt)))

# Set up the initial wavefunction and potential
hydrogen_string = make_hydrogen_string(num_gauss_potential, external_field_params)
try:
    initial_linear_coeffs, params_init = read_best_coefficients(
        num_gauss_wavefunction=num_gauss_wavefunction,
        num_gauss_potential=num_gauss_potential,
    )
except FileNotFoundError:
    if dt.imag:
        print(
            f"No ground state data for ng_wf={ngwf}, ng_pot={ngp}. "
            f"Using random initial state for imaginary time propagation."
        )
        D = 3  # Hydrogen is 3D
        widths = np.logspace(-2, 2, num_gauss_wavefunction)
        params_init = jnp.zeros((num_gauss_wavefunction, D, 4), dtype=jnp.float64)
        params_init = params_init.at[:, :, 0].set(widths[:, None])
        initial_linear_coeffs = jnp.ones(num_gauss_wavefunction, dtype=jnp.complex128)
    else:
        raise FileNotFoundError(
            f"No ground state data for ng_wf={ngwf}, ng_pot={ngp}. Cannot proceed with real time propagation without initial state."
        )
hydrogen_wf = generalPotentialSolver(params_init, hydrogen_string)

# Set up the Rothe error and gradient functions
rothe_error, rothe_gradient, SHH2 = set_up_functions(
    params_init,
    hydrogen_string,
    splitting_type="none",  # We ignore splitting type not none for now
)

# Set up Rothe solver
rothesolver = RotheSolver(
    SHH2,
    dt,
    0,
    params_init,
    initial_linear_coeffs,
    rothe_grad_jax=rothe_gradient,
    rothe_nograd=rothe_error,
    splitting_type="none",
    output_config=OutputConfig(
        name=string_name,
        polynomial_string=hydrogen_string,
        out_dir="./wave_function_data",
    ),
)
nsteps = resume_or_initialize_rothe(
    rothesolver,
    nsteps,
    rothesolver.name,
    "none",
    hydrogen_string,
    dt,
    initial_error=0,
    params=params_init,
    coeffs=initial_linear_coeffs,
    t_start=t_start,
)

rothesolver.propagate(nsteps, maxiter=maxiter)
