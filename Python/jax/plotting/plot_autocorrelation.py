#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fftfreq, ifft

sys.path.insert(0, str(Path(__file__).parent.parent))

from rothe.io import load_rothe_file


def parse_sim_name(raw: str) -> str:
    """Normalize CLI input to a bare simulation name.

    Accepts either:
      - bare sim names (e.g. "Foo__epsilon=1e-3__lambda=1e-6")
      - full filenames (e.g. "wave_function_data/Foo__...__none.h5")
    """
    # If a path is provided, only keep the basename.
    name = Path(raw).name

    # Strip .h5 extension when user passes a full filename.
    if name.endswith(".h5"):
        name = name[:-3]

    # If a trailing splitting token is present ("__none", "__strang", ...),
    # drop it so load_rothe_file(sim_name, "none") resolves correctly.
    head, sep, tail = name.rpartition("__")
    if sep and "=" not in tail:
        return head

    return name


def load(filename: str) -> tuple[np.ndarray, np.ndarray]:
    """Load ``t`` and ``double_autocorrelation`` from a simulation file.

    Input:
        - ``filename``: full ``.h5`` path or simulation name-like token

    Output:
        - ``times``: NumPy array of time points
        - ``autocorr``: NumPy complex array of double autocorrelation values

    If autocorrelation is non-finite at ``t=0``, it is set to ``1+0j``.
    """
    sim_name = parse_sim_name(filename)
    data = load_rothe_file(sim_name, "none")

    times = np.asarray(data.get("t", []))
    autocorr_list = data["double_autocorrelation"]

    autocorr = np.array([v if v is not None else np.nan + 0j for v in autocorr_list], dtype=np.complex128)

    finite_mask = np.isfinite(np.real(autocorr)) & np.isfinite(np.imag(autocorr))
    t0_mask = np.isclose(np.abs(times), 0.0)
    fix_mask = (~finite_mask) & t0_mask
    autocorr[fix_mask] = 1.0 + 0.0j

    return times, autocorr


def sort_by_time(times: np.ndarray, autocorr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sort ``times`` and reorder ``autocorr`` with the same argsort index."""
    idx = np.argsort(times)
    return times[idx], autocorr[idx]


def calculate_spectrum(times: np.ndarray, autocorr: np.ndarray):
    # Use legacy convention explicitly: x = 2t in damping.
    x = 2.0 * np.asarray(times)
    y = np.array(autocorr) * np.exp(-x / 30)
    yf = np.asarray(ifft(y))
    N = len(times)
    # If t_final=100, this gives T_final=200 as requested.
    T = (2.0 * times[-1]) / N
    xf = fftfreq(N, T)[: N // 2] * np.pi
    yf = np.real(yf[: N // 2])
    return xf, yf


def main():
    parser = argparse.ArgumentParser(description="Read time points and double autocorrelation from simulation files")
    parser.add_argument(
        "sim_names", nargs="+", help="One or more simulation names (or full filenames in wave_function_data)"
    )
    args = parser.parse_args()
    all_times = []
    all_autocorrs = []
    for raw_name in args.sim_names:
        sim_name = parse_sim_name(raw_name)
        times, autocorr = load(raw_name)
        times, autocorr = sort_by_time(times, autocorr)
        all_times.append(times)
        all_autocorrs.append(autocorr)
        print(f"Loaded {sim_name!r}:")
    for times, autocorr in zip(all_times, all_autocorrs):
        # xf, yf = calculate_spectrum(times, autocorr)
        plt.plot(2 * times, abs(autocorr) ** 2, label="Real part")
    plt.show()


if __name__ == "__main__":
    main()
