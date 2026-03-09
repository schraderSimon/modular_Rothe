"""
Extract ground states from imaginary-time propagation data.

Scans wave_function_data/ for imaginary-time HDF5 files, computes the energy
variance at each time step, picks the minimum, and writes the best
coefficients to Hydrogen_ground_state_data/{name}.npz.

Usage:  python scripts/extract_groundstates.py
"""

import pickle
import re
from pathlib import Path

import numpy as np

import jax.numpy as jnp
from rothe.io import load_rothe_file
from rothe.systems.common import ExternalFieldParams
from rothe.systems.hydrogen import make_hydrogen_string
from rothe.wavefunction import generalPotentialSolver

# ── paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CACHE_FILE = PROJECT_ROOT / "results" / "variance_cache.pkl"
DATA_DIRS = [
    PROJECT_ROOT / "wave_function_data",
    PROJECT_ROOT / "wave_function_data" / "imaginary_time_propagation_data_old",
]
OUTPUT_PATH = PROJECT_ROOT / "Hydrogen_ground_state_data"


# ── helpers ────────────────────────────────────────────────────────────────
def parse_simulation_name(sim_name):
    """Return (n_wf, n_pot) from either old or new naming format."""
    m = re.match(r"Hydrogen_(\d+)_(\d+)", sim_name)
    if m:
        return int(m.group(1)), int(m.group(2))
    wf = re.search(r"ng_wf=(\d+)", sim_name)
    pot = re.search(r"ng_pot=(\d+)", sim_name)
    if wf and pot:
        return int(wf.group(1)), int(pot.group(1))
    raise ValueError(f"Cannot parse: {sim_name}")


def is_imaginary_time_simulation(sim_name):
    if re.match(r"Hydrogen_\d+_\d+$", sim_name):
        return True
    return bool(re.search(r"dt=[^_]*i", sim_name))


def _variance_and_energy(params, coeffs, poly_string):
    """<H^2> - <H>^2 and <H> for one time step."""
    S, H, H2 = generalPotentialSolver(params, params, poly_string).calculate_SHH2(t=0, splitting_type="none")
    c = coeffs / jnp.sqrt(jnp.vdot(coeffs, S @ coeffs))
    E = jnp.vdot(c, H @ c)
    return float(jnp.real(jnp.vdot(c, H2 @ c) - E**2)), float(jnp.real(E))


# ── variance cache ─────────────────────────────────────────────────────────
def load_variance_cache():
    if CACHE_FILE.exists():
        with open(CACHE_FILE, "rb") as f:
            return pickle.load(f)
    return {}


def save_variance_cache(cache):
    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE_FILE, "wb") as f:
        pickle.dump(cache, f)


def _compute_variances_for_file(sim_name, directory, max_time=-500j):
    """Return (times, variances, energies) for one simulation file."""
    n_wf, n_pot = parse_simulation_name(sim_name)
    data = load_rothe_file(sim_name, "none", path=str(directory))
    poly = make_hydrogen_string(n_pot, ExternalFieldParams())

    times = data["t"]
    # Use absolute imaginary depth so the filter works for both +dt*j and -dt*j.
    mask = np.abs(np.imag(times)) >= abs(np.imag(max_time))
    idx = np.where(mask)[0]
    times = times[mask]
    var = np.empty(len(idx))
    eng = np.empty(len(idx))

    print(f"  {sim_name}: {len(idx)} steps ...", flush=True)
    for i, j in enumerate(idx):
        var[i], eng[i] = _variance_and_energy(data["params"][j], data["coeffs"][j], poly)
    return times, var, eng


def compute_all_variances(max_time=100j):
    """Compute (or load cached) variances for every imaginary-time simulation."""
    cache = load_variance_cache()
    modified = False
    for d in DATA_DIRS:
        if not d.exists():
            continue
        for f in sorted(d.iterdir()):
            if not f.name.endswith("__none.h5"):
                continue
            name = f.name.replace("__none.h5", "")
            if not is_imaginary_time_simulation(name) or name in cache:
                continue
            try:
                cache[name] = _compute_variances_for_file(name, d, max_time)
                modified = True
            except Exception as e:
                print(f"  {name}: ERROR - {e}")
    if modified:
        save_variance_cache(cache)
    return cache


# ── extract & save best coefficients ───────────────────────────────────────
def extract_best_coefficients(cache=None):
    if cache is None:
        cache = compute_all_variances()
    OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

    for sim_name, (times, variances, energies) in sorted(cache.items()):
        try:
            n_wf, n_pot = parse_simulation_name(sim_name)
        except ValueError:
            continue

        best = int(np.argmin(variances))

        data = None
        for d in DATA_DIRS:
            try:
                data = load_rothe_file(sim_name, "none", path=str(d))
                break
            except FileNotFoundError:
                pass
        if data is None:
            continue

        j = int(np.argmin(np.abs(data["t"] - times[best])))
        if np.abs(data["t"][j] - times[best]) > 1e-10:
            continue

        np.savez(
            OUTPUT_PATH / f"Hydrogen_{n_wf}_{n_pot}.npz",
            linear_coeffs=data["coeffs"][j],
            params=data["params"][j],
            variance=variances[best],
            energy=energies[best],
            time=times[best],
            num_gauss_wavefunction=n_wf,
            num_gauss_potential=n_pot,
        )
        print(f"  Hydrogen_{n_wf}_{n_pot}: var={variances[best]:.2e}, E={energies[best]:.6f}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract ground states from imaginary-time propagation data.")
    parser.add_argument(
        "file",
        nargs="?",
        default=None,
        help="Path to a specific HDF5 file to process (e.g. wave_function_data/mysim__none.h5). "
        "If omitted, all imaginary-time files in DATA_DIRS are scanned.",
    )
    args = parser.parse_args()

    if args.file is not None:
        path = Path(args.file).resolve()
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        sim_name = path.name.replace("__none.h5", "")
        # Use max_time=0j so all steps are included regardless of how far the run got.
        cache = {sim_name: _compute_variances_for_file(sim_name, path.parent, max_time=0j)}
    else:
        cache = compute_all_variances()

    extract_best_coefficients(cache)
