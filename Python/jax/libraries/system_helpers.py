"""
Hydrogen and Henon-Heiles system helpers
"""

from pathlib import Path

import numpy as np
import pandas as pd
from attr import dataclass

from jax import numpy as jnp
from libraries.general_wf import generalPotentialSolver
from libraries.jax_Rothe import setUpRotheErrorAndGradient_jit


@dataclass
class ExternalFieldParams:
    E0: float = 0.06  # field strength in atomic units
    omega: float = 0.057  # frequency in atomic units
    N_T: int = 3  # Number of cycles

    @property
    def intensity(self) -> float:
        return 3.5094e16 * self.E0**2

    @property
    def t_final(self) -> float:
        return 2 * jnp.pi / self.omega * self.N_T


def set_up_functions(params_init, potential_string, splitting_type="none"):
    """Set up Rothe error/gradient functions and SHH2 function.

    Returns:
    - rothe_error: error function
    - rothe_vg_jit: gradient function
    - SHH2: function to compute S, H, H2 matrices
    """
    rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type)
    system = generalPotentialSolver(params_init, potential_string)

    def SHH2(t, params, splitting_type=splitting_type):
        system.update_parameters(params)
        return system.calculate_SHH2(t=t, splitting_type=splitting_type)

    return rothe_error, rothe_vg_jit, SHH2


def eval_initial_state_properties(initial_linear_coeffs, params_init, SHH2):
    initial_linear_coeffs = jnp.array(initial_linear_coeffs)
    S, H, H2 = SHH2(0, params_init)
    norm = jnp.vdot(initial_linear_coeffs, S @ initial_linear_coeffs)
    initial_linear_coeffs = initial_linear_coeffs / jnp.sqrt(norm)
    E0 = jnp.vdot(initial_linear_coeffs, H @ initial_linear_coeffs)
    Esq = jnp.vdot(initial_linear_coeffs, H2 @ initial_linear_coeffs)
    print("Initial energy:", E0)
    print("Initial norm:", norm)
    print("Initial variance:", Esq - E0**2)


def make_hydrogen_string(
    num_gauss_potential: int, external_field_params: ExternalFieldParams
) -> str:
    """
    Construct the potential string for the hydrogen atom
    for a given number of Gaussian functions in the potential expansion.
    """
    _project_root = Path(__file__).parent.parent
    csv_path = (
        _project_root / "gaussian_Coulomb" / f"coeffs_mu=100_N={num_gauss_potential}.csv"
    )
    df = pd.read_csv(csv_path)
    lincoeffs = df["linear"]
    width = df["nonlinear"]
    hydrogen_string = """
    dimension 3
    polynomial
    x2: E_0*sin(omega*t)*sin(omega*t/N_T/2)**2
    exponential
    """
    E0 = external_field_params.E0
    omega = external_field_params.omega
    N_T = external_field_params.N_T
    for lc, w in zip(lincoeffs, width):
        hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"  # Centered at the origin
    hydrogen_string += (
        f"constants\nN_T: {N_T}\nomega: {omega}\nE_0: {E0}"  # Time dependent field
    )
    return hydrogen_string


def set_up_hydrogen_string_name(
    external_field_params: ExternalFieldParams, ngp: int, ngwf: int, dt: complex
) -> str:
    E0 = external_field_params.E0
    omega = external_field_params.omega
    N_T = external_field_params.N_T

    if dt.real and dt.imag:
        raise ValueError("dt must be either purely real or purely imaginary.")

    val, suffix = (dt.real, "") if dt.real else (dt.imag, "i")
    dt = float(val)
    dt_string = f"{dt:.2f}{suffix}"
    if dt.imag:
        assert (
            E0 == 0.0
        ), "For imaginary time propagation, the external field strength E0 must be zero."
    string_name = f"hydrogen__E0={E0}__omega={omega}__N_T={N_T}__ng_pot={ngp}__ng_wf={ngwf}_dt={dt_string}"
    return string_name


def initialize_henon_heiles_params(n: int, D: int, seed: int = 42):
    """
    Initialize parameters for Henon-Heiles system.

    Parameters:
    - n: number of Gaussians
    - D: dimension
    - seed: random seed

    Returns:
    - params_init: initialized parameters array of shape (n, D, 4)
    """
    np.random.seed(seed)
    params_init = jnp.zeros((n, D, 4), dtype=jnp.float64)
    params_init = params_init.at[:, :, 0].set(1 / jnp.sqrt(2))  # Width parameters
    params_init = params_init.at[0, :, 2].set(2.0)  # Set mu to (2,...,2)
    for i in range(1, n):
        params_init = params_init.at[i, :, 2].set(
            2 + np.random.uniform(-0.5, 0.5, (D,))
        )  # Move to the right
    return params_init


def initialize_henon_heiles_coeffs(n: int, params_init, potential_string):
    """
    Initialize coefficients for Henon-Heiles system.

    Parameters:
    - n: number of Gaussians
    - params_init: initialized parameters
    - potential_string: the potential string

    Returns:
    - coeffs_init: normalized coefficients
    """
    oscillator = generalPotentialSolver(params_init, potential_string)
    Smat = oscillator.calculate_S()
    coeffs_init = jnp.array([1.0] * n, dtype=jnp.complex128)
    coeffs_init = jnp.zeros_like(coeffs_init)
    coeffs_init = coeffs_init.at[0].set(1.0)
    norm = coeffs_init.T @ Smat @ coeffs_init
    coeffs_init = coeffs_init / jnp.sqrt(norm)

    return coeffs_init


def resume_or_initialize_rothe(
    rothesolver,
    nsteps_total,
    solver_name,
    splitting_type,
    potential_string,
    dt,
    initial_error=0.0,
    params=None,
    coeffs=None,
    t_start=None,
):
    """
    Resume Rothe solver from checkpoint or initialize from scratch.

    Parameters:
    - rothesolver: RotheSolver instance
    - nsteps_total: total number of steps to simulate
    - solver_name: name of the solver (for file identification)
    - splitting_type: splitting type used
    - potential_string: the potential string
    - dt: time step
    - initial_error: initial error value (default 0.0)
    - params: initial parameters
    - coeffs: initial coefficients
    - t_start: if given, resume from this time (or the latest checkpoint <= t_start)

    Returns:
    - nsteps_remaining: number of steps left to propagate
    """
    from libraries.file_handling import save_rothe_state

    try:
        info = rothesolver.resume_from_file(t_start=t_start)
        print(f"Resumed at t={info['t']} (step {info['idx']}), trimmed={info['trimmed']}")
        return nsteps_total - info["idx"]
    except FileNotFoundError:
        save_rothe_state(
            solver_name,
            splitting_type,
            potential_string,
            dt,
            0.0,
            float(initial_error),
            params=params,
            coeffs=coeffs,
        )
        print("No prior file. Starting fresh at t=0.")
        return nsteps_total
