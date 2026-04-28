#!/usr/bin/env python3
import argparse
import re

# Add parent directory to path for imports
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).parent.parent))

from rothe.io import load_rothe_file

OLD_DATA_DIR = Path(__file__).parent.parent / "old_data"


def parse_sim_input(raw: str) -> tuple[str, str, Path | None]:
    """Return (sim_name, splitting_type, path_override) from a CLI input token.

    Accepts any of:
    - plain sim name: ``hydrogen__...__num_dynamic=10``
    - sim name with explicit splitting token: ``hydrogen__...__num_dynamic=10__none``
    - full/relative HDF5 path: ``wave_function_data/hydrogen__...__num_dynamic=10__none.h5``
    """
    token_path = Path(raw)
    token = raw
    path_override = None

    if token_path.suffix == ".h5":
        token = token_path.stem
        if token_path.parent != Path("."):
            path_override = token_path.parent
    elif token_path.parent != Path("."):
        token = token_path.name
        path_override = token_path.parent

    parts = token.split("__")
    if len(parts) >= 2 and "=" not in parts[-1]:
        splitting_type = parts[-1]
        sim_name = "__".join(parts[:-1])
    else:
        splitting_type = "none"
        sim_name = token

    return sim_name, splitting_type, path_override


def parse_E0_from_sim_name(sim_name: str) -> float | None:
    """Extract E0 from a new-format sim name, e.g. 'hydrogen__E0=0.06__...'."""
    m = re.search(r"E0=([0-9.]+)", sim_name)
    return float(m.group(1)) if m else None


def _extract_param(sim_name: str, key: str) -> str | None:
    """Extract a single key=value token from a double-underscore simulation name."""
    m = re.search(rf"(?:^|__){re.escape(key)}=([^_]+)(?:__|$)", sim_name)
    return m.group(1) if m else None


def make_plot_label(sim_name: str) -> str:
    """Compact label: only Gaussian count and epsilon."""
    n_gauss = _extract_param(sim_name, "ng_wf")
    if n_gauss is None:
        n_gauss = _extract_param(sim_name, "num_dynamic")
    eps = _extract_param(sim_name, "epsilon")

    pieces = []
    if n_gauss is not None:
        pieces.append(f"N={n_gauss}")
    if eps is not None:
        pieces.append(f"ε={eps}")

    return ", ".join(pieces) if pieces else sim_name


def find_ref_file(E0: float, old_data_dir: Path) -> tuple[Path, float] | None:
    """Return (path, epsilon) for the old_data file matching E0 with the lowest epsilon."""
    candidates = []
    for p in sorted(old_data_dir.glob("*.h5")):
        m_e0 = re.search(r"_E0([0-9.]+)_", p.name)
        if m_e0 is None:
            continue
        if not np.isclose(float(m_e0.group(1)), E0, rtol=1e-3):
            continue
        m_eps = re.search(r"_epsilon([0-9.e+\-]+)_", p.name)
        if m_eps is None:
            continue
        candidates.append((float(m_eps.group(1)), p))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0])
    return candidates[0][1], candidates[0][0]


def load_old_ref(filepath: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load (times, Re(dpm)) from an old-format HDF5 file."""
    with h5py.File(filepath, "r") as f:
        times = np.array(f["times"])
        dpm = np.real(np.array(f["dpm"]))
    return times, dpm


def main():
    parser = argparse.ArgumentParser(description="Plot z-dipole moment evolution")
    parser.add_argument(
        "sim_names", nargs="+", help="One or more simulation names (or full filenames in wave_function_data)"
    )
    parser.add_argument("--output", default=None, help="Output file path (default: ./figures/dipole.png)")
    parser.add_argument(
        "--omega0",
        type=float,
        default=0.057,
        help="Fundamental frequency ω₀ in a.u. (default: 0.057, corresponding to 800 nm)",
    )
    parser.add_argument(
        "--no-ref",
        dest="plot_ref",
        action="store_false",
        default=True,
        help="Disable plotting the lowest-epsilon reference from old_data for matching E0",
    )
    args = parser.parse_args()

    # Ensure figures directory exists
    figures_dir = Path(__file__).parent.parent / "figures"
    figures_dir.mkdir(exist_ok=True)

    # Determine output path
    if args.output is None:
        output_path = figures_dir / "dipole.png"
    else:
        output_path = Path(args.output)
        if output_path.parent == Path("."):
            output_path = figures_dir / output_path.name

    fig, (ax_dip, ax_hhg) = plt.subplots(1, 2, figsize=(14, 5))

    time_label = "Real Time (t)"
    plotted_ref_E0s: set[float] = set()

    for raw_name in args.sim_names:
        sim_name, splitting_type, path_override = parse_sim_input(raw_name)
        plot_label = make_plot_label(sim_name)

        try:
            data = load_rothe_file(sim_name, splitting_type, path=str(path_override) if path_override else None)
        except FileNotFoundError:
            print(
                f"Error: Could not find data for input {raw_name!r} "
                f"(resolved to sim_name={sim_name!r}, splitting={splitting_type!r}) — skipping."
            )
            continue

        times = data["t"]
        dipole_list = data.get("dipole", None)

        if dipole_list is None or all(v is None for v in dipole_list):
            print(f"Error: No dipole moment data found in {sim_name!r} — skipping.")
            continue

        # Build array from list, filling None entries with NaN
        D_dip = next(v.shape[-1] for v in dipole_list if v is not None)
        dipole = np.full((len(dipole_list), D_dip), np.nan)
        for i, v in enumerate(dipole_list):
            if v is not None:
                dipole[i] = v

        # Use only the last (z) component
        dz = np.real(dipole[:, -1])

        # Determine time axis
        if np.iscomplexobj(times):
            times_plot = np.imag(times)
            time_label = "Imaginary Time (Im(t))"
        else:
            times_plot = np.real(times)
            time_label = "Real Time (t)"

        # Guard against partially-written files where t and dipole have different lengths
        n = min(len(times_plot), len(dz))
        if n < max(len(times_plot), len(dz)):
            print(
                f"  Warning: t ({len(times_plot)}) and dipole ({len(dz)}) lengths differ "
                f"for {sim_name!r}; truncating to {n}."
            )
        times_plot = times_plot[:n]
        dz = dz[:n]

        ax_dip.plot(times_plot, dz, label=plot_label)

        # HHG spectrum — drop NaN entries (initial t=0 state has no dipole)
        valid = np.isfinite(dz)
        if valid.sum() < 2:
            print(f"  Warning: not enough finite dipole samples for HHG in {sim_name!r} — skipping spectrum.")
            continue
        omega, spec = compute_hhg_spectrum(times_plot[valid], dz[valid])
        # Keep only positive frequencies; rescale to harmonic order n = ω/ω₀
        pos = omega > 0
        harmonic_order = omega[pos] / args.omega0
        ax_hhg.semilogy(harmonic_order, omega[pos] ** 2 * spec[pos], label=plot_label)

        # Optionally overlay the lowest-epsilon reference from old_data
        if args.plot_ref and OLD_DATA_DIR.exists():
            E0 = parse_E0_from_sim_name(sim_name)
            if E0 is not None and E0 not in plotted_ref_E0s:
                result = find_ref_file(E0, OLD_DATA_DIR)
                if result is None:
                    print(f"  No reference found in old_data for E0={E0}.")
                else:
                    ref_path, ref_eps = result
                    ref_label = f"ref E0={E0} (ε={ref_eps:.2e})"
                    ref_times, ref_dz = load_old_ref(ref_path)
                    ax_dip.plot(ref_times[:n], ref_dz[:n], linestyle="--", color="gray", alpha=0.7, label=ref_label)
                    ref_omega, ref_spec = compute_hhg_spectrum(ref_times[:n], ref_dz[:n])
                    ref_pos = ref_omega > 0
                    ref_ho = ref_omega[ref_pos] / args.omega0
                    ax_hhg.semilogy(
                        ref_ho,
                        ref_omega[ref_pos] ** 2 * ref_spec[ref_pos],
                        linestyle="--",
                        color="gray",
                        alpha=0.7,
                        label=ref_label,
                    )
                    plotted_ref_E0s.add(E0)
                    print(f"  Reference: {ref_path.name}")

    ax_dip.set_xlabel(time_label)
    ax_dip.set_ylabel("Dipole moment $d_z$")
    ax_dip.set_title("z-component of dipole moment")
    if ax_dip.lines:
        ax_dip.legend()
    ax_dip.grid(True, alpha=0.3)

    ax_hhg.set_xlabel(r"Harmonic order $n = \omega/\omega_0$")
    ax_hhg.set_ylabel(r"$\omega^2 |d(\omega)|^2$ (arb. units)")
    ax_hhg.set_title(f"HHG spectrum (velocity form, $\\omega_0={args.omega0}$ a.u.)")
    ax_hhg.set_xlim(1, 50)
    if ax_hhg.lines:
        ax_hhg.legend()
    ax_hhg.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved to {output_path}")
    plt.show()


def compute_hhg_spectrum(time_points, dipole_moment, hann_window=True):
    """Compute the HHG spectrum from a dipole moment time series."""

    dip_np = np.array(dipole_moment)
    dip_detrended = scipy.signal.detrend(dip_np, type="constant")

    if hann_window:
        window = np.sin(np.pi * np.array(time_points) / time_points[-1]) ** 2
        dip_windowed = dip_detrended * window
        spec = np.abs(scipy.fftpack.fftshift(scipy.fftpack.fft(dip_windowed))) ** 2
    else:
        spec = np.abs(scipy.fftpack.fftshift(scipy.fftpack.fft(dip_detrended))) ** 2

    dt = time_points[1] - time_points[0]
    freq = scipy.fftpack.fftshift(scipy.fftpack.fftfreq(len(time_points)))
    omega_arr = freq * 2 * np.pi / dt
    return omega_arr, spec


if __name__ == "__main__":
    main()
