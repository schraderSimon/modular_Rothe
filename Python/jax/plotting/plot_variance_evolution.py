"""
Plot variance evolution as a function of imaginary time for Hydrogen ground state calculations.

Reads pre-computed variance data from the cache produced by
``initial_state_calculations/groundstate_from_imaginary_data.py``
and creates one subplot per wave-function basis size.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from initial_state_calculations.groundstate_from_imaginary_data import (
    compute_all_variances,
    load_variance_cache,
    parse_simulation_name,
)

COLORS = [
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


def main():
    # Try loading existing cache; compute only if empty
    cache = load_variance_cache()
    if not cache:
        print("No cached variance data found — computing now …")
        cache = compute_all_variances()

    # Group cache entries by n_wf → [(n_pot, sim_name, times, variances, energies)]
    groups = {}
    for sim_name, (times, variances, energies) in cache.items():
        try:
            n_wf, n_pot = parse_simulation_name(sim_name)
        except ValueError:
            continue
        groups.setdefault(n_wf, []).append((n_pot, sim_name, times, variances, energies))

    for n_wf in groups:
        groups[n_wf].sort(key=lambda x: x[0])  # sort by n_pot

    sorted_n_wf = sorted(groups.keys())
    n_panels = len(sorted_n_wf)

    if n_panels == 0:
        print("Nothing to plot (cache is empty).")
        return

    fig, axes = plt.subplots(
        1, n_panels, figsize=(n_panels * 7, 8), sharex=True, squeeze=False
    )
    axes = axes[0]  # squeeze=False gives 2-D array

    for ax, n_wf in zip(axes, sorted_n_wf):
        for i, (n_pot, sim_name, times, variances, _energies) in enumerate(groups[n_wf]):
            imag_times = np.imag(times)
            if len(imag_times) == 0:
                continue
            ax.semilogy(
                -imag_times,
                variances,
                label=f"{n_pot} Gaussians for potential",
                color=COLORS[i % len(COLORS)],
                linewidth=1.5,
            )

        ax.set_ylabel("Variance", fontsize=12)
        ax.set_ylim(1e-11, 1e-1)
        ax.set_title(f"Wave function: {n_wf} Gaussians", fontsize=14)
        ax.legend(loc="upper right", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 80)
        ax.set_xlabel("Imaginary time (−Im(t))", fontsize=12)

    plt.suptitle(
        "Energy Variance in Imaginary Time\n(Hydrogen Ground State)",
        fontsize=14,
        y=0.98,
    )
    plt.tight_layout()

    figures_dir = Path(__file__).parent.parent / "figures"
    figures_dir.mkdir(exist_ok=True)
    output_file = figures_dir / "variance_evolution.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {output_file}")
    plt.show()


if __name__ == "__main__":
    main()
