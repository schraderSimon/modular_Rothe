#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare grid, Python-Rothe, and (optionally) Fortran-Rothe solutions
===================================================================

Put this script in your Python/ directory and run it from there.
"""

from pathlib import Path
import numpy as np
import h5py
import matplotlib.pyplot as plt


# ─────────────────────────── utilities ────────────────────────────
def evaluate_basis(p, c, x):
    """Reconstruct ψ(x) from a Gaussian basis (row-wise p = [a, b, μ, p])."""
    psi = np.zeros_like(x, dtype=complex)
    for (a, b, mu, p_param), ci in zip(p, c):
        A = a**2 + 1j * b
        N = (2 * A.real / np.pi) ** 0.25
        psi += ci * N * np.exp(-A * (x - mu) ** 2 + 1j * p_param * (x - mu))
    return psi


def nearest_index(t_ref, t_query):
    """Return indices in t_ref that are closest to each element of t_query."""
    return np.abs(t_ref[:, None] - t_query).argmin(axis=0)


# ────────────────────────── grid reference ─────────────────────────
grid_data = np.load("output/wavefunction.npz")
x       = grid_data["x"]
t_grid  = grid_data["t"]
psi_ref = grid_data["psi"]          # shape (N_t, N_x)
dx      = x[1] - x[0]

# ─────────────────────── Python-Rothe trajectory ───────────────────
with h5py.File("HO_1D_trajectory.h5", "r") as f:
    keys = sorted(k for k in f if k.startswith("step_"))
    n_frames = len(keys)
    psi_py   = np.empty((n_frames, x.size), dtype=complex)
    err_py   = np.empty(n_frames)

    for i, key in enumerate(keys):
        grp = f[key]
        p   = grp["p"][...]
        c   = grp["c"][...]
        if i==10:
            print(p)
            print(c)
        psi = evaluate_basis(p, c, x)
        psi_py[i] = psi
        err_py[i] = np.sqrt(np.sum(np.abs(psi - psi_ref[i])**2) * dx)

# ───────────────────── optional Fortran trajectory ─────────────────
ft_dir  = Path(__file__).parent.parent / "Fortran" / "output"
p_file  = ft_dir / "nonlinear_data.npy"
c_file  = ft_dir / "coefficients.npy"
t_file  = ft_dir / "times.npy"

has_fortran = all(p.exists() for p in (p_file, c_file, t_file))

if has_fortran:
    p_out = np.load(p_file)  # shape (n_g, N_b, n_tF)
    print(p_out.shape)
    print(p_out[:,:,10])
    c_out = np.load(c_file)        # shape (n_g,     N_tF)
    print(c_out)
    print(c_out.shape)
    t_f   = np.load(t_file)        # (N_tF,)
    print(t_f)
    n_tF  = t_f.shape[-1]
    print(n_tF)
    psi_ft = np.empty((n_tF, x.size), dtype=complex)
    err_ft = np.empty(n_tF)

    # Map each Fortran time-point to the nearest grid time
    idx_map = nearest_index(t_grid, t_f)

    for j in range(1,n_tF):
        psi_j = evaluate_basis(p_out[:, :, j-1], c_out[:, j-1], x)
        psi_ft[j] = psi_j
        # L² error versus *grid* reference at the matching time
        g_idx  = idx_map[j]
        err_ft[j] = np.sqrt(np.sum(np.abs(psi_j - psi_ref[g_idx])**2) * dx)
    print(err_ft)
# ──────────────────────────── PLOTTING ────────────────────────────
plt.figure()
plt.plot(t_grid, err_py, label="Python Rothe", linewidth=2)

if has_fortran:
    plt.plot(t_f, err_ft, ":", label="Fortran Rothe", linewidth=2)
    print(
        f"Fortran data found in {ft_dir}. "
        "Plotted true L² error vs. grid reference."
    )
else:
    print("No Fortran output found – skipping second curve.")

plt.xlabel("Time")
plt.ylabel("L² error to grid")
plt.title("Gaussian-basis L² Error")
plt.legend()

# ─────────────────────── final-time comparison ────────────────────
plt.figure()
plt.plot(x, psi_ref[-2].real, label="Grid")
plt.plot(x, psi_py[-2].real, "--", label="Basis (Python)")
if has_fortran:
    plt.plot(x, psi_ft[-1].real, ":", label="Basis (Fortran)")
plt.xlabel("x")
plt.ylabel(r"Re $\psi(x,\,t_{\mathrm{final}})$")
plt.title("Real part at final time")
plt.legend()
plt.tight_layout()
plt.show()
