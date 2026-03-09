#!/usr/bin/env python3
import argparse

# Add parent directory to path for imports
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from rothe.io import load_rothe_file


def parse_sim_name(raw: str) -> str:
    """Strip a file extension token appended by double-underscore, if present."""
    parts = raw.split("__")
    if parts and "=" not in parts[-1]:
        return "__".join(parts[:-1])
    return raw


def main():
    parser = argparse.ArgumentParser(description="Plot absolute double autocorrelation evolution")
    parser.add_argument(
        "sim_names", nargs="+", help="One or more simulation names (or full filenames in wave_function_data)"
    )
    parser.add_argument("--output", default=None, help="Output file path (default: ./figures/autocorrelation.png)")
    args = parser.parse_args()

    # Ensure figures directory exists
    figures_dir = Path(__file__).parent.parent / "figures"
    figures_dir.mkdir(exist_ok=True)

    # Determine output path
    if args.output is None:
        output_path = figures_dir / "autocorrelation.png"
    else:
        output_path = Path(args.output)
        if output_path.parent == Path("."):
            output_path = figures_dir / output_path.name

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    time_label = "Real Time (t)"

    for raw_name in args.sim_names:
        sim_name = parse_sim_name(raw_name)

        try:
            data = load_rothe_file(sim_name, "none")
        except FileNotFoundError:
            print(f"Error: Could not find data for {sim_name!r} — skipping.")
            continue

        times = data["t"]
        autocorr_list = data.get("double_autocorrelation", None)

        if autocorr_list is None or all(v is None for v in autocorr_list):
            print(f"Error: No double autocorrelation data found in {sim_name!r} — skipping.")
            continue

        autocorr = np.array([v if v is not None else np.nan + 0j for v in autocorr_list])
        abs_autocorr = np.abs(autocorr)

        # Determine time axis
        if np.iscomplexobj(times):
            times_plot = np.imag(times)
            time_label = "Imaginary Time (Im(t))"
        else:
            times_plot = np.real(times)
            time_label = "Real Time (t)"

        # Guard against partially-written files where t and autocorrelation have different lengths
        n = min(len(times_plot), len(abs_autocorr))
        if n < max(len(times_plot), len(abs_autocorr)):
            print(
                f"  Warning: t ({len(times_plot)}) and double_autocorrelation ({len(abs_autocorr)}) lengths differ "
                f"for {sim_name!r}; truncating to {n}."
            )
        times_plot = times_plot[:n]
        abs_autocorr = abs_autocorr[:n]

        ax.plot(times_plot, abs_autocorr, label=sim_name)

    ax.set_xlabel(time_label)
    ax.set_ylabel(r"$|C(2t)| = |\langle \Psi^*(t)|\Psi(t)\rangle|$")
    ax.set_title("Absolute double autocorrelation")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
