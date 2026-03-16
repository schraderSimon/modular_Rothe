"""
Hydrogen atom system helpers.

Builds the potential string for a Coulomb potential expanded in Gaussians,
optionally with a time-dependent external field.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from rothe._config import PROJECT_ROOT

from .common import ExternalFieldParams

# Default data directory for pre-computed ground state coefficients
_GROUND_STATE_DIR = PROJECT_ROOT / "Hydrogen_ground_state_data"


def make_hydrogen_string(num_gauss_potential: int, external_field_params: ExternalFieldParams) -> str:
    """
    Construct the potential string for the hydrogen atom
    for a given number of Gaussian functions in the potential expansion.
    """
    csv_path = PROJECT_ROOT / "gaussian_Coulomb" / f"coeffs_mu=100_N={num_gauss_potential}.csv"
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
        hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"
    hydrogen_string += f"constants\nN_T: {N_T}\nomega: {omega}\nE_0: {E0}"
    return hydrogen_string


def make_hydrogen_squared_string(num_gauss_potential: int, num_gauss_fit: int) -> str:
    """
    Build a potential string for the pre-fitted approximation of V_gauss^2.

    Reads ``square_coeffs_mu=100_Norig={num_gauss_potential}_Nfit={num_gauss_fit}.csv``
    and returns an exponential-only string that can be passed as ``squared_string``
    to ``generalPotentialSolver``.
    """
    csv_path = (
        PROJECT_ROOT / "gaussian_Coulomb" / f"square_coeffs_mu=100_Norig={num_gauss_potential}_Nfit={num_gauss_fit}.csv"
    )
    df = pd.read_csv(csv_path)
    lincoeffs = df["linear"]
    width = df["nonlinear"]
    D = 3
    string = "dimension 3\nexponential\n"
    for lc, w in zip(lincoeffs, width):
        string += f"{lc}, {w}, {[0.0] * D}\n"
    return string


def set_up_hydrogen_string_name(
    external_field_params: ExternalFieldParams,
    ngp: int,
    ngwf: int,
    dt: complex,
    epsilon: float,
    regularization_lambda: float,
    num_dynamic: int,
) -> str:
    """Build a descriptive filename stem for a hydrogen simulation run."""
    E0 = external_field_params.E0
    omega = external_field_params.omega
    N_T = external_field_params.N_T

    if dt.real and dt.imag:
        raise ValueError("dt must be either purely real or purely imaginary.")

    val, suffix = (dt.real, "") if dt.real else (dt.imag, "i")
    dt_float = float(val)
    dt_string = f"{dt_float:.2f}{suffix}"
    if dt.imag:
        assert E0 == 0.0, "For imaginary time propagation, E0 must be zero."
    string_name = (
        f"hydrogen__E0={E0}__omega={omega}__N_T={N_T}__ng_pot={ngp}__ng_wf={ngwf}"
        f"_dt={dt_string}__epsilon={float(epsilon):.3e}__lambda={float(regularization_lambda):.3e}"
        f"__num_dynamic={int(num_dynamic)}"
    )
    return string_name


def read_best_coefficients(num_gauss_potential, num_gauss_wavefunction, data_dir=None):
    """Load (linear_coeffs, params) from a previously extracted .npz file.

    Parameters:
        num_gauss_potential: number of Gaussians in the potential expansion
        num_gauss_wavefunction: number of Gaussians in the wavefunction
        data_dir: directory containing the .npz files (default: Hydrogen_ground_state_data/)
    """
    if data_dir is None:
        data_dir = _GROUND_STATE_DIR
    data_dir = Path(data_dir)
    path = data_dir / f"Hydrogen_{num_gauss_wavefunction}_{num_gauss_potential}.npz"
    if not path.exists():
        avail = (
            sorted((int(f.stem.split("_")[1]), int(f.stem.split("_")[2])) for f in data_dir.glob("Hydrogen_*.npz"))
            if data_dir.exists()
            else []
        )
        raise FileNotFoundError(f"{path.name} not found. Available: {avail}")
    data = np.load(path)
    return data["linear_coeffs"], data["params"]
