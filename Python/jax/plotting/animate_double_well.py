#!/usr/bin/env python3

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


def parse_args():
    parser = argparse.ArgumentParser(description="Animate density evolution in the 1D double well")
    parser.add_argument("input_file", help="Path to the saved HDF5 trajectory")
    parser.add_argument("--points", type=int, default=1000, help="Number of spatial grid points")
    parser.add_argument("--interval", type=int, default=40, help="Base delay between frames in ms")
    parser.add_argument("--speed", type=float, default=2.0, help="Playback speed multiplier")
    parser.add_argument("--repeat", action="store_true", help="Loop the animation")
    parser.add_argument("--xmin", type=float, default=-4.0, help="Left edge of the spatial window")
    parser.add_argument("--xmax", type=float, default=4.0, help="Right edge of the spatial window")
    return parser.parse_args()


def load_trajectory(filepath: Path):
    with h5py.File(filepath, "r") as handle:
        polynomial_string = handle["polynomial_string"][()]
        if isinstance(polynomial_string, bytes):
            polynomial_string = polynomial_string.decode("utf-8")

        steps_group = handle["steps"]
        step_keys = sorted(steps_group.keys())
        times = []
        params = []
        coeffs = []

        for key in step_keys:
            step = steps_group[key]
            times.append(step.attrs["t"])
            params.append(step["params"][...])
            coeffs.append(step["coeffs"][...])

    return polynomial_string, np.asarray(times), params, coeffs


def evaluate_density_1d(x_grid: np.ndarray, params_t: np.ndarray, coeffs_t: np.ndarray) -> np.ndarray:
    if params_t.ndim != 3 or params_t.shape[1] != 1:
        raise ValueError("This animator only supports 1D trajectories.")

    a = params_t[:, 0, 0]
    b = params_t[:, 0, 1]
    mu = params_t[:, 0, 2]
    p = params_t[:, 0, 3]

    dx = x_grid[None, :] - mu[:, None]
    normalization = (2.0 * a**2 / np.pi) ** 0.25
    basis_values = normalization[:, None] * np.exp(-(a**2 + 1j * b)[:, None] * dx**2 + 1j * p[:, None] * dx)
    psi = coeffs_t @ basis_values
    return np.abs(psi) ** 2


def build_spatial_grid(x_min: float, x_max: float, num_points: int) -> np.ndarray:
    return np.linspace(x_min, x_max, num_points)


def double_well_potential(x_grid: np.ndarray) -> np.ndarray:
    return 0.25 * x_grid**4 - 0.5 * x_grid**2


def main():
    args = parse_args()
    input_path = Path(args.input_file).expanduser().resolve()

    if args.xmin >= args.xmax:
        raise ValueError("Expected xmin < xmax.")
    if args.speed <= 0:
        raise ValueError("Expected speed > 0.")

    if not input_path.exists():
        raise FileNotFoundError(f"Could not find input file: {input_path}")

    polynomial_string, times, params_history, coeffs_history = load_trajectory(input_path)
    if len(times) == 0:
        raise ValueError(f"No saved steps found in {input_path}")

    x_grid = build_spatial_grid(args.xmin, args.xmax, args.points)
    density_history = [
        evaluate_density_1d(x_grid, params_t, coeffs_t) for params_t, coeffs_t in zip(params_history, coeffs_history)
    ]
    density_max = max(float(np.max(density)) for density in density_history)
    frame_interval = max(1, int(round(args.interval / args.speed)))

    fig, ax_density = plt.subplots(figsize=(10, 5))
    (density_line,) = ax_density.plot([], [], color="tab:blue", lw=2)
    ax_density.set_xlim(float(x_grid[0]), float(x_grid[-1]))
    ax_density.set_ylim(0.0, 1.05 * density_max if density_max > 0 else 1.0)
    ax_density.set_xlabel("x")
    ax_density.set_ylabel(r"Density $|\psi(x,t)|^2$")
    ax_density.grid(True, alpha=0.3)

    ax_potential = ax_density.twinx()
    potential = double_well_potential(x_grid)
    ax_potential.plot(x_grid, potential, color="0.4", linestyle="--", alpha=0.8)
    ax_potential.set_ylabel("V(x)", color="0.4")
    ax_potential.tick_params(axis="y", colors="0.4")

    if np.iscomplexobj(times):
        time_values = np.imag(times)
        time_label = "Im(t)"
    else:
        time_values = np.real(times)
        time_label = "t"

    title_text = ax_density.set_title("")
    if "x0x0x0x0: 0.25" not in polynomial_string or "x0x0: -0.5" not in polynomial_string:
        title_text.set_text("Density evolution")

    def init():
        density_line.set_data([], [])
        return density_line, title_text

    def update(frame_index: int):
        density_line.set_data(x_grid, density_history[frame_index])
        title_text.set_text(f"Density evolution, {time_label} = {time_values[frame_index]:.6g}")
        return density_line, title_text

    animation = FuncAnimation(
        fig,
        update,
        frames=len(density_history),
        init_func=init,
        interval=frame_interval,
        blit=True,
        repeat=args.repeat,
        cache_frame_data=False,
    )

    fig.tight_layout()
    plt.show()
    return animation


if __name__ == "__main__":
    main()
