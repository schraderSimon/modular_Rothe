"""
Plot the dipole moment evolution from a Hydrogen propagation simulation.

Usage:
    python3 plot_dipole_moment.py --num_gauss_potential 11 --num_gauss_wavefunction 15
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from libraries.file_handling import load_rothe_file


def plot_dipole_moment(num_gauss_potential, num_gauss_wavefunction, output_file=None):
    """
    Load and plot the dipole moment evolution for a Hydrogen simulation.

    Args:
        num_gauss_potential: Number of Gaussians for the potential
        num_gauss_wavefunction: Number of Gaussians for the wavefunction
        output_file: Optional file to save the plot (e.g., 'dipole_moment.png')
    """
    sim_name = f"Hydrogen_{num_gauss_wavefunction}_{num_gauss_potential}"

    try:
        data = load_rothe_file(sim_name, "none", path="./wave_function_data")
    except FileNotFoundError:
        print(f"Error: Could not find data for {sim_name}")
        print("Make sure the simulation has been run and the HDF5 file exists.")
        return

    times = data["t"]
    dipole = data.get("dipole")

    if dipole is None:
        print(f"Error: No dipole moment data found in {sim_name}__none.h5")
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
        axes = [axes]  # Make it iterable

    # Plot each dimension
    for d in range(D):
        ax = axes[d]
        ax.plot(times_plot, dipole[:, d], "o-", linewidth=2, markersize=4)
        ax.set_xlabel(time_label, fontsize=12)
        ax.set_ylabel(f"⟨{dim_labels[d]}⟩", fontsize=12)
        ax.set_title(f"Dipole Moment: {dim_labels[d]}-component", fontsize=12)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"Dipole Moment Evolution: Hydrogen with {num_gauss_wavefunction} WF Gaussians, "
        f"{num_gauss_potential} Potential Gaussians",
        fontsize=14,
        y=1.00,
    )
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
        print(f"Figure saved to {output_file}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot dipole moment evolution from Rothe simulation"
    )
    parser.add_argument(
        "--num_gauss_potential",
        type=int,
        required=True,
        help="Number of Gaussians for the potential",
    )
    parser.add_argument(
        "--num_gauss_wavefunction",
        type=int,
        required=True,
        help="Number of Gaussians for the wavefunction",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional: save figure to this file instead of showing it",
    )

    args = parser.parse_args()

    plot_dipole_moment(
        args.num_gauss_potential, args.num_gauss_wavefunction, output_file=args.output
    )


if __name__ == "__main__":
    main()
