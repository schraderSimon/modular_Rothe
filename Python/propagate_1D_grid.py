#!/usr/bin/env python3
"""
Brute-force Crank–Nicolson propagation of the 1D Schrödinger equation
for a shifted Gaussian in a harmonic oscillator (ω=1).
Uses dense matrices only (no sparse or Thomas solver).
"""
import numpy as np
from pathlib import Path

def main():
    # Grid and time parameters
    x_min, x_max = -10.0, 10.0
    dx = 0.02
    x = np.arange(x_min, x_max + dx, dx)
    N = x.size

    dt = 0.01
    t_max = 10.0
    nt = int(t_max / dt)
    times = np.linspace(0.0, t_max, nt + 1)

    # Initial normalized Gaussian centered at x=1
    a=1/np.sqrt(3)
    b=0.1
    mu=1
    q= 0
    psi1 = np.exp(-(a**2+1j*b)*(x-mu)**2+1j*q*(x-mu))
    psi1 /= np.sqrt(np.trapezoid(np.abs(psi1)**2, x))
    psi2= np.exp(-(a**2+1j*b)*(x+mu)**2+1j*q*(x+mu))
    psi2 /= np.sqrt(np.trapezoid(np.abs(psi2)**2, x))
    psi3 = np.exp(-(a**2+1j*b)*(x-2*mu)**2+1j*q*(x-2*mu))
    psi3 /= np.sqrt(np.trapezoid(np.abs(psi3)**2, x))
    psi4 = np.exp(-(a**2+1j*b)*(x+2*mu)**2+1j*q*(x+2*mu))
    psi4 /= np.sqrt(np.trapezoid(np.abs(psi4)**2, x))
    psi=psi1+psi2+psi3+psi4
    # Construct kinetic operator T (dense)
    T = np.zeros((N, N), dtype=complex)
    off = -1.0 / (2 * dx**2)
    diag = 1.0 / (dx**2)
    for i in range(N):
        T[i, i] = diag
        if i > 0:
            T[i, i-1] = off
        if i < N - 1:
            T[i, i+1] = off

    # Potential operator V
    V = 0.5 * x**2
    Vmat = np.diag(V)

    # Total Hamiltonian H = T + V
    H = T + Vmat
    # Dirichlet BC: enforce ψ=0 at boundaries by zeroing H rows/cols
    H[0, :] = 0; H[:, 0] = 0; H[0, 0] = 0
    H[-1, :] = 0; H[:, -1] = 0; H[-1, -1] = 0

    # Crank-Nicolson matrices
    I = np.eye(N, dtype=complex)
    A = I + 1j * dt / 2 * H
    B = I - 1j * dt / 2 * H
    # Enforce BC rows in A and B
    for idx in (0, N-1):
        A[idx, :] = 0; A[idx, idx] = 1
        B[idx, :] = 0; B[idx, idx] = 1

    # Allocate storage
    psi_store = np.empty((nt + 1, N), dtype=complex)
    psi_store[0] = psi

    # Time evolution
    for n in range(nt):
        d = B @ psi
        psi = np.linalg.solve(A, d)
        psi[0] = psi[-1] = 0.0
        psi_store[n + 1] = psi

    # Save results
    out_dir = Path('output')
    out_dir.mkdir(exist_ok=True)
    np.savez(out_dir / 'wavefunction.npz', x=x, t=times, psi=psi_store)
    print(f"Saved wavefunctions to {out_dir/'wavefunctions.npz'}")

if __name__ == '__main__':
    main()
