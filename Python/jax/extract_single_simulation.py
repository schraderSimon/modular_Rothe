"""
Extract best coefficients from a specific Hydrogen simulation H5 file
and save to a .npz file in Hydrogen_ground_state_data directory.
"""

import os

import numpy as np

from libraries.file_handling import load_rothe_file

# Paths
DATA_PATH = "./wave_function_data"
OUTPUT_PATH = "./Hydrogen_ground_state_data"


def extract_specific_simulation(sim_name):
    """
    Extract the best coefficients from a specific simulation and save to .npz file.

    Args:
        sim_name: Name of the simulation (e.g., "Hydrogen_50_21")
    """
    print(f"Processing {sim_name}...")

    # Parse simulation name
    parts = sim_name.split("_")
    if len(parts) < 3 or parts[0] != "Hydrogen":
        raise ValueError(f"Invalid simulation name format: {sim_name}")

    n_wf = int(parts[1])  # number of Gaussians in wave function
    n_pot = int(parts[2])  # number of Gaussians in potential

    # Load the wave function data
    try:
        data = load_rothe_file(sim_name, "none", path=DATA_PATH)
        print(f"  Loaded data with {len(data['t'])} time points")
    except FileNotFoundError as e:
        print(f"  ERROR: wave function data not found: {e}")
        return False

    # Calculate variances for all time points
    print("  Calculating variances...")
    times = data["t"]
    energies = data["energy"]
    variances = data["variance"]

    # Find minimum variance index
    min_idx = np.argmin(variances)
    min_variance = variances[min_idx]
    min_time = times[min_idx]
    min_energy = energies[min_idx]

    print(
        f"  Best state at t={min_time:.3f}: variance={min_variance:.2e}, energy={min_energy:.6f}"
    )

    # Extract coefficients and parameters for the best state
    linear_coeffs = data["coeffs"][min_idx]
    params = data["params"][min_idx]

    print(f"  Linear coeffs shape: {linear_coeffs.shape}")
    print(f"  Params shape: {params.shape}")

    # Create output directory
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    # Save to .npz file
    output_file = os.path.join(OUTPUT_PATH, f"{sim_name}.npz")
    np.savez(
        output_file,
        linear_coeffs=linear_coeffs,
        params=params,
        variance=min_variance,
        energy=min_energy,
        time=min_time,
        num_gauss_wavefunction=n_wf,
        num_gauss_potential=n_pot,
    )

    print(f"  SUCCESS: Saved to {output_file}")
    return True


if __name__ == "__main__":
    # Extract the specific simulation you want
    sim_name = "Hydrogen_50_21"

    print("=" * 60)
    print(f"Extracting ground state for {sim_name}")
    print("=" * 60)

    success = extract_specific_simulation(sim_name)

    if success:
        print("\nDone! You can now use read_best_coefficients.py to load the data.")
        print(
            f"Example: read_best_coefficients(num_gauss_potential=21, num_gauss_wavefunction=50)"
        )
    else:
        print("\nFailed to extract ground state data.")
