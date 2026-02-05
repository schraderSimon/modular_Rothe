"""
Extract the best (lowest variance) coefficients from all Hydrogen simulations
and save them to individual .npz files in Hydrogen_ground_state_data directory.
"""

import os
import pickle

import numpy as np

from libraries.file_handling import load_rothe_file

# Paths
CACHE_FILE = "results/variance_cache.pkl"
DATA_PATH = "./wave_function_data"
OUTPUT_PATH = "./Hydrogen_ground_state_data"


def extract_and_save_all():
    """
    For each simulation in the variance cache, extract the best coefficients
    and save them to a .npz file.
    """
    # Load variance cache
    if not os.path.exists(CACHE_FILE):
        raise FileNotFoundError(f"Variance cache file not found: {CACHE_FILE}")

    with open(CACHE_FILE, "rb") as f:
        cache = pickle.load(f)

    # Create output directory
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    print(f"Extracting best coefficients from {len(cache)} simulations...")
    print("-" * 70)

    for sim_name, (times, variances, energies) in cache.items():
        if not sim_name.startswith("Hydrogen_"):
            continue

        # Parse simulation name
        parts = sim_name.split("_")
        n_wf = int(parts[1])
        n_pot = int(parts[2])

        # Find minimum variance index
        min_idx = np.argmin(variances)
        min_variance = variances[min_idx]
        min_time = times[min_idx]
        min_energy = energies[min_idx]

        # Load the wave function data
        try:
            data = load_rothe_file(sim_name, "none", path=DATA_PATH)
        except FileNotFoundError:
            print(f"  {sim_name}: wave function data not found, skipping")
            continue

        # Find matching time in full data
        all_times = data["t"]
        time_diffs = np.abs(all_times - min_time)
        data_idx = np.argmin(time_diffs)

        if time_diffs[data_idx] > 1e-10:
            print(f"  {sim_name}: time mismatch, skipping")
            continue

        # Extract coefficients and parameters
        linear_coeffs = data["coeffs"][data_idx]
        params = data["params"][data_idx]

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

        print(
            f"  {sim_name}: variance={min_variance:.2e}, energy={min_energy:.6f} -> {output_file}"
        )

    print("-" * 70)
    print(f"Done! Files saved to {OUTPUT_PATH}/")


if __name__ == "__main__":
    extract_and_save_all()
