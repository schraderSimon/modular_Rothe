#!/usr/bin/env python3
import argparse

# Add parent directory to path for imports
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from libraries.file_handling import load_rothe_file


def main():
    parser = argparse.ArgumentParser(description="Plot dipole moment evolution")
    parser.add_argument(
        "sim_name", help="Simulation name (or full filename in wave_function_data)"
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output file path (default: ./figures/{sim_name}_dipole.png)",
    )
    args = parser.parse_args()

    # Smart filename parsing: if last part after split doesn't contain "=", it's an extension
    # e.g., "henon_heiles__mu=100__N=10.h5" -> "henon_heiles__mu=100__N=10"
    parts = args.sim_name.split("__")
    if parts and "=" not in parts[-1]:
        # Remove the extension part
        sim_name = "__".join(parts[:-1])
    else:
        sim_name = args.sim_name

    # Ensure figures directory exists
    figures_dir = Path(__file__).parent.parent / "figures"
    figures_dir.mkdir(exist_ok=True)

    # Determine output path
    if args.output is None:
        output_path = figures_dir / f"{sim_name}_dipole.png"
    else:
        output_path = Path(args.output)
        if output_path.parent == Path("."):
            output_path = figures_dir / output_path.name

    # Read data
    try:
        data = load_rothe_file(sim_name, "none", path="../wave_function_data")
    except FileNotFoundError:
        print(f"Error: Could not find data for {sim_name}")
        print("Make sure the simulation has been run and the HDF5 file exists.")
        return

    # Extract time steps and dipole moments
    times = data["t"]
    dipole = data.get("dipole", None)

    # Check if dipole data exists
    if dipole is None:
        print(f"Error: No dipole moment data found in {sim_name}")
        print(
            "Make sure the simulation was run with the updated code that stores dipole moments."
        )
        return

    # Extract dimension labels
    D = dipole.shape[1]
    dim_labels = ["x", "y", "z"][:D]

    # Determine if time is real or imaginary
    is_imaginary = np.iscomplexobj(times)
    if is_imaginary:
        times_plot = np.imag(times)
        time_label = "Imaginary Time (Im(t))"
    else:
        times_plot = np.real(times)
        time_label = "Real Time (t)"

    # Create figure
    fig, axes = plt.subplots(1, D, figsize=(5 * D, 5), sharex=True)
    if D == 1:
        axes = np.array([axes])

    # Plot each dimension
    for i, (ax, label) in enumerate(zip(axes, dim_labels)):
        ax.plot(times_plot, np.real(dipole[:, i]), label="Real part")
        ax.plot(times_plot, np.imag(dipole[:, i]), label="Imag part")
        ax.set_ylabel("Dipole moment")
        ax.set_title(f"{label}-component")
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.xlabel(time_label)
    plt.tight_layout()

    # Save figure
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
