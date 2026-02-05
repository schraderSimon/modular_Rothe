"""
Plot variance evolution as a function of imaginary time for Hydrogen ground state calculations.

Creates 3 subplots (one per number of Gaussians: 15, 25, 30) showing how the variance
evolves over imaginary time propagation for different potential Gaussian basis sizes.
"""

import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import jax.numpy as jnp
from libraries.file_handling import load_rothe_file
from libraries.general_wf import generalPotentialSolver

# Cache file for storing computed variances
CACHE_FILE = "results/variance_cache.pkl"


def calculate_variance_at_timestep(params, coeffs, polynomial_string):
    """
    Calculate the variance <H²> - <H>² for a given set of parameters and coefficients.

    Args:
        params: Array of shape (n, D, 4) with Gaussian parameters
        coeffs: Array of shape (n,) with complex linear coefficients
        polynomial_string: String defining the potential

    Returns:
        variance: The energy variance
        energy: The expectation value <H>
    """
    solver = generalPotentialSolver(params, polynomial_string)
    S, H, H2 = solver.calculate_SHH2(t=0, splitting_type="none")

    # Normalize coefficients
    norm = jnp.vdot(coeffs, S @ coeffs)
    coeffs_normalized = coeffs / jnp.sqrt(norm)

    # Calculate expectation values
    energy = jnp.vdot(coeffs_normalized, H @ coeffs_normalized)
    energy_sq = jnp.vdot(coeffs_normalized, H2 @ coeffs_normalized)

    variance = jnp.real(energy_sq - energy**2)
    return float(variance), float(jnp.real(energy))


def get_hydrogen_potential_string(num_gauss_potential):
    """
    Load the Gaussian Coulomb potential approximation from file.

    Args:
        num_gauss_potential: Number of Gaussians used to approximate 1/r

    Returns:
        polynomial_string: String defining the potential
    """
    df = pd.read_csv(f"gaussian_Coulomb/coeffs_mu=100_N={num_gauss_potential}.csv")
    lincoeffs = df["linear"]
    width = df["nonlinear"]

    hydrogen_string = """
    dimension 3
    exponential
    """
    for lc, w in zip(lincoeffs, width):
        hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"

    return hydrogen_string


def load_variance_cache():
    """
    Load cached variance data from file if it exists.

    Returns:
        dict: Cache dictionary mapping sim_name to (times, variances, energies)
    """
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, "rb") as f:
            return pickle.load(f)
    return {}


def save_variance_cache(cache):
    """
    Save variance cache to file.

    Args:
        cache: Dictionary mapping sim_name to (times, variances, energies)
    """
    with open(CACHE_FILE, "wb") as f:
        pickle.dump(cache, f)


def process_simulation_file(
    sim_name, splitting_type, path="./wave_function_data", max_time=-100j, cache=None
):
    """
    Load simulation data and calculate variance at each time step.
    Uses caching to avoid recalculating if data already exists.

    Args:
        sim_name: Name of the simulation (e.g., "Hydrogen_25_31")
        splitting_type: The splitting type used (e.g., "none")
        path: Path to the wave function data directory
        max_time: Maximum imaginary time to process (e.g., -100j)
        cache: Dictionary for caching results (optional)

    Returns:
        times: Array of imaginary times
        variances: Array of variances at each time
        energies: Array of energies at each time
    """
    # Check cache first
    if cache is not None and sim_name in cache:
        print(f"Loading {sim_name} from cache...")
        return cache[sim_name]

    # Parse simulation name to get number of Gaussians for potential
    parts = sim_name.split("_")
    num_gauss_wf = int(parts[1])
    num_gauss_pot = int(parts[2])

    # Load data
    data = load_rothe_file(sim_name, splitting_type, path=path)

    # Get polynomial string from the potential Gaussian basis
    polynomial_string = get_hydrogen_potential_string(num_gauss_pot)

    # Arrays to store results
    all_times = data["t"]
    params_all = data["params"]
    coeffs_all = data["coeffs"]

    # Filter to only process times up to max_time
    imag_times = np.imag(all_times)
    max_imag = np.imag(max_time)  # This is negative (e.g., -100)
    mask = imag_times >= max_imag  # Keep times with imag part >= -100

    indices = np.where(mask)[0]
    n_steps = len(indices)

    times = all_times[mask]
    variances = np.zeros(n_steps)
    energies = np.zeros(n_steps)

    print(f"Processing {sim_name}: {n_steps} time steps (up to t={max_time})...")

    for i, idx in enumerate(indices):
        var, energy = calculate_variance_at_timestep(
            params_all[idx], coeffs_all[idx], polynomial_string
        )
        variances[i] = var
        energies[i] = energy

        if (i + 1) % 20 == 0 or i == n_steps - 1:
            print(f"  Step {i+1}/{n_steps}, t={times[i]:.2f}, var={var:.6f}, E={energy:.6f}")

    result = (times, variances, energies)

    # Store in cache
    if cache is not None:
        cache[sim_name] = result

    return result


def main():
    """
    Create a figure with 3 subplots showing variance evolution for different
    numbers of Gaussians in the wave function (15, 25, 30).
    """
    # Configuration
    num_gauss_wf_list = [15, 25, 30, 35, 50]
    max_time = -100j  # Only consider up to this imaginary time

    # Find all available files
    data_path = "./wave_function_data"
    files = os.listdir(data_path)

    # Group files by number of wave function Gaussians
    file_groups = {n: [] for n in num_gauss_wf_list}

    for f in files:
        if f.endswith(".h5") and f.startswith("Hydrogen_"):
            parts = f.replace("__none.h5", "").split("_")
            n_wf = int(parts[1])
            n_pot = int(parts[2])
            if n_wf in num_gauss_wf_list:
                file_groups[n_wf].append((n_pot, f"Hydrogen_{n_wf}_{n_pot}"))

    # Sort each group by number of potential Gaussians
    for n in num_gauss_wf_list:
        file_groups[n].sort(key=lambda x: x[0])

    # Load variance cache
    cache = load_variance_cache()
    cache_modified = False

    # Create figure with 3 subplots
    fig, axes = plt.subplots(
        1, len(num_gauss_wf_list), figsize=(len(num_gauss_wf_list) * 7, 8), sharex=True
    )

    # Color map for different potential Gaussian numbers
    # colors = plt.cm.viridis(np.linspace(0.1, 0.9, 10))
    colors = [
        "blue",
        "red",
        "green",
        "orange",
        "purple",
        "brown",
        "pink",
        "gray",
        "olive",
        "cyan",
    ]
    for ax_idx, (ax, n_wf) in enumerate(zip(axes, num_gauss_wf_list)):
        print(f"\n=== Processing wave function with {n_wf} Gaussians ===")

        for color_idx, (n_pot, sim_name) in enumerate(file_groups[n_wf]):
            try:
                print("Hmm")
                was_cached = sim_name in cache
                times, variances, energies = process_simulation_file(
                    sim_name, "none", path=data_path, max_time=max_time, cache=cache
                )
                if not was_cached:
                    cache_modified = True
                print("Hemm")
                # Times are already filtered in process_simulation_file
                imag_times = np.imag(times)

                if len(imag_times) == 0:
                    print(f"  No data points within time range for {sim_name}")
                    continue
                print("homm")
                # Plot variance vs negative imaginary time
                ax.semilogy(
                    -imag_times,
                    variances,
                    label=f"{n_pot} Gaussians to approximate erf(100r)/r",
                    color=colors[color_idx % len(colors)],
                    linewidth=1.5,
                )

            except Exception as e:
                print(f"  Error processing {sim_name}: {e}")

        ax.set_ylabel("Variance", fontsize=12)
        ax.set_ylim(1e-11, 1e-1)
        ax.set_title(f"Wave function with {n_wf} Gaussians", fontsize=14)
        ax.legend(loc="upper right", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 80)
        ax.set_xlabel("Imaginary time (-Im(t))", fontsize=12)
    # Set common x-label

    plt.suptitle(
        "Energy Variance Evolution in Imaginary Time\n(Hydrogen Ground State)",
        fontsize=14,
        y=0.98,
    )
    plt.tight_layout()

    # Save cache if modified
    if cache_modified:
        save_variance_cache(cache)
        print(f"\nVariance cache saved to {CACHE_FILE}")

    # Save figure
    output_file = "variance_evolution.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {output_file}")

    plt.show()


if __name__ == "__main__":
    main()
