"""
Utility to read the best (lowest variance) linear and nonlinear coefficients
for a given Hydrogen simulation from pre-extracted .npz files.
"""

import os

import numpy as np

# Default path to pre-extracted ground state data
GROUND_STATE_PATH = "./Hydrogen_ground_state_data"


def read_best_coefficients(
    num_gauss_potential,
    num_gauss_wavefunction,
    data_path=GROUND_STATE_PATH,
):
    """
    Load the linear and nonlinear coefficients with the lowest variance
    for a given Hydrogen simulation from pre-extracted .npz files.

    Args:
        num_gauss_potential: Number of Gaussians used to approximate the Coulomb potential
        num_gauss_wavefunction: Number of Gaussians in the wave function basis
        data_path: Path to the ground state data directory

    Returns:
        linear_init: Array of shape (num_gauss_wavefunction,) with complex linear coefficients
        params_init: Array of shape (num_gauss_wavefunction, 3, 4) with Gaussian parameters

    Raises:
        FileNotFoundError: If the ground state file is not found
    """
    sim_name = f"Hydrogen_{num_gauss_wavefunction}_{num_gauss_potential}"
    filepath = os.path.join(data_path, f"{sim_name}.npz")

    if not os.path.exists(filepath):
        available = list_available_simulations(data_path)
        available_str = [f"({n_wf}, {n_pot})" for n_wf, n_pot, _ in available]
        raise FileNotFoundError(
            f"Ground state file not found: {filepath}. "
            f"Available (n_wf, n_pot): {available_str}"
        )

    data = np.load(filepath)
    linear_init = data["linear_coeffs"]
    params_init = data["params"]

    return linear_init, params_init


def get_best_variance_info(
    num_gauss_potential,
    num_gauss_wavefunction,
    data_path=GROUND_STATE_PATH,
):
    """
    Get detailed information about the best (lowest variance) state.

    Args:
        num_gauss_potential: Number of Gaussians used to approximate the Coulomb potential
        num_gauss_wavefunction: Number of Gaussians in the wave function basis
        data_path: Path to the ground state data directory

    Returns:
        dict with keys:
            - linear_coeffs: Complex linear coefficients
            - params: Gaussian parameters
            - variance: The minimum variance value
            - energy: The energy at the minimum variance state
            - time: The imaginary time at which minimum variance occurs
            - sim_name: The simulation name
    """
    sim_name = f"Hydrogen_{num_gauss_wavefunction}_{num_gauss_potential}"
    filepath = os.path.join(data_path, f"{sim_name}.npz")

    if not os.path.exists(filepath):
        available = list_available_simulations(data_path)
        available_str = [f"({n_wf}, {n_pot})" for n_wf, n_pot, _ in available]
        raise FileNotFoundError(
            f"Ground state file not found: {filepath}. "
            f"Available (n_wf, n_pot): {available_str}"
        )

    data = np.load(filepath)

    return {
        "linear_coeffs": data["linear_coeffs"],
        "params": data["params"],
        "variance": float(data["variance"]),
        "energy": float(data["energy"]),
        "time": complex(data["time"]),
        "sim_name": sim_name,
    }


def list_available_simulations(data_path=GROUND_STATE_PATH):
    """
    List all available Hydrogen simulations in the ground state data directory.

    Returns:
        list of tuples: (num_gauss_wavefunction, num_gauss_potential, min_variance)
    """
    if not os.path.exists(data_path):
        return []

    results = []
    for filename in os.listdir(data_path):
        if filename.startswith("Hydrogen_") and filename.endswith(".npz"):
            # Parse filename: Hydrogen_{n_wf}_{n_pot}.npz
            name = filename.replace(".npz", "")
            parts = name.split("_")
            n_wf = int(parts[1])
            n_pot = int(parts[2])

            # Load variance from file
            filepath = os.path.join(data_path, filename)
            data = np.load(filepath)
            min_var = float(data["variance"])

            results.append((n_wf, n_pot, min_var))

    # Sort by num_gauss_wavefunction, then num_gauss_potential
    results.sort(key=lambda x: (x[0], x[1]))
    return results


if __name__ == "__main__":
    # Example usage
    print("Available simulations:")
    print("-" * 60)
    print(f"{'n_wf':>6} {'n_pot':>6} {'min_variance':>15}")
    print("-" * 60)

    for n_wf, n_pot, min_var in list_available_simulations():
        print(f"{n_wf:>6} {n_pot:>6} {min_var:>15.2e}")

    print("\n" + "=" * 60)
    print("Example: reading best coefficients for n_pot=11, n_wf=15")
    print("=" * 60)

    try:
        linear_init, params_init = read_best_coefficients(
            num_gauss_potential=11, num_gauss_wavefunction=15
        )
        print(f"linear_init shape: {linear_init.shape}")
        print(f"params_init shape: {params_init.shape}")

        info = get_best_variance_info(num_gauss_potential=11, num_gauss_wavefunction=15)
        print(f"Minimum variance: {info['variance']:.6e}")
        print(f"Energy at min variance: {info['energy']:.6f}")
        print(f"Time at min variance: {info['time']}")
    except (FileNotFoundError, KeyError) as e:
        print(f"Error: {e}")
