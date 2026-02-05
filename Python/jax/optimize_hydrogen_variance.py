import numpy as np
import pandas as pd
from scipy.optimize import minimize

import jax
import jax.numpy as jnp
from libraries.general_wf import generalPotentialSolver
from libraries.utils import calculate_variance

# Build the hydrogen potential string from stored Gaussian fit coefficients
# The linear coefficients are negated to represent the attractive Coulomb tail.
df = pd.read_csv("gaussian_Coulomb/coeffs_mu=100_N=21.csv")
lincoeffs = df["linear"]
widths = df["nonlinear"]

HYDROGEN_STRING = """
    dimension 3
    exponential
    """
for lc, w in zip(lincoeffs, widths):
    HYDROGEN_STRING += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"


# =============================================================================
# Core helper functions (shared by all optimization routines)
# =============================================================================


def params_from_log_widths(log_w: jnp.ndarray, n: int) -> jnp.ndarray:
    """Convert log-widths to the params array expected by generalPotentialSolver."""
    w = jnp.exp(jnp.asarray(log_w))
    params = jnp.zeros((n, 3, 4), dtype=jnp.float64)
    return params.at[:, :, 0].set(w[:, None])


def compute_orthogonalized_matrices(log_w: jnp.ndarray, n: int, compute_H2: bool = False):
    """Compute S, H (and optionally H2) matrices and orthogonalize them.

    Returns:
        H_orth: Orthogonalized Hamiltonian
        H2_orth: Orthogonalized H² (or None if compute_H2=False)
    """
    params = params_from_log_widths(log_w, n)
    system = generalPotentialSolver(params, HYDROGEN_STRING)
    S = system.calculate_S()
    H = system.calculate_H()
    H2 = system.calculate_H2() if compute_H2 else None

    # Regularize and symmetrize overlap matrix
    S = 0.5 * (S + S.conj().T) + 1e-8 * jnp.eye(S.shape[0])

    # Orthogonalize via Löwdin transformation
    w, U = jnp.linalg.eigh(S)
    X = U * (1.0 / jnp.sqrt(w))[None, :]
    H_orth = X.conj().T @ H @ X
    H2_orth = X.conj().T @ H2 @ X if compute_H2 else None

    return H_orth, H2_orth


def compute_ground_state(H_orth: jnp.ndarray):
    """Diagonalize H_orth and return ground state energy and eigenvector."""
    eigvals, eigvecs = jnp.linalg.eigh(H_orth)
    energy = jnp.real(eigvals[0])
    gs = eigvecs[:, 0]
    return energy, gs


def compute_variance_and_energy(log_w: jnp.ndarray, n: int) -> tuple[float, float]:
    """Compute variance and energy for a given set of log-widths."""
    H_orth, H2_orth = compute_orthogonalized_matrices(log_w, n, compute_H2=True)
    energy, gs = compute_ground_state(H_orth)
    variance = float(calculate_variance(gs, H2_orth, float(energy)))
    return variance, float(energy)


# =============================================================================
# Objective functions for optimization
# =============================================================================


def make_energy_objective(n: int):
    """Create a JIT-compiled energy objective function and its gradient."""

    def energy_only(log_w):
        H_orth, _ = compute_orthogonalized_matrices(log_w, n, compute_H2=False)
        energy, _ = compute_ground_state(H_orth)
        return energy

    return jax.jit(jax.value_and_grad(energy_only))


def make_variance_objective(n: int):
    """Create a JIT-compiled variance objective function and its gradient."""

    def variance_only(log_w):
        H_orth, H2_orth = compute_orthogonalized_matrices(log_w, n, compute_H2=True)
        energy, gs = compute_ground_state(H_orth)
        expectation_H2 = jnp.real(gs.conj() @ H2_orth @ gs)
        return expectation_H2 - energy**2

    return jax.jit(jax.value_and_grad(variance_only))


def make_logspace_energy_objective(n: int):
    """Create energy objective over (minval, maxval) parameters for log-spaced widths."""

    def energy_from_bounds(bounds):
        minval, maxval = bounds
        log_w = jnp.log(jnp.logspace(minval, maxval, n))
        H_orth, _ = compute_orthogonalized_matrices(log_w, n, compute_H2=False)
        energy, _ = compute_ground_state(H_orth)
        return energy

    return jax.jit(jax.value_and_grad(energy_from_bounds))


def make_logspace_variance_objective(n: int):
    """Create variance objective over (minval, maxval) parameters for log-spaced widths."""

    def variance_from_bounds(bounds):
        minval, maxval = bounds
        log_w = jnp.log(jnp.logspace(minval, maxval, n))
        H_orth, H2_orth = compute_orthogonalized_matrices(log_w, n, compute_H2=True)
        energy, gs = compute_ground_state(H_orth)
        expectation_H2 = jnp.real(gs.conj() @ H2_orth @ gs)
        return expectation_H2 - energy**2

    def log_variance_from_bounds(bounds):
        return jnp.log(variance_from_bounds(bounds))

    return jax.jit(jax.value_and_grad(log_variance_from_bounds))


# =============================================================================
# Main optimization routines
# =============================================================================


def optimize_for_n(
    n: int,
    log_width_init: np.ndarray | None = None,
    objective: str = "variance",
    maxiter: int = 0,
) -> dict:
    """Optimize Gaussian widths for the hydrogen ground state.

    Args:
        n: Basis size (number of Gaussians)
        log_width_init: Initial log-widths. If None, use log-spaced grid from 0.1 to 10.
        objective: What to minimize - "energy" or "variance"
        maxiter: Maximum optimization iterations

    Returns:
        Dictionary with optimization results including variance, energy, and optimal widths.
    """
    if log_width_init is None:
        log_width_init = np.log(np.logspace(-1.0, 1, n))

    # Select objective function
    if objective == "energy":
        objective_vg = make_energy_objective(n)
    elif objective == "variance":
        objective_vg = make_variance_objective(n)
    else:
        raise ValueError(f"Unknown objective: {objective}. Use 'energy' or 'variance'.")

    def scipy_objective(log_w_np):
        val, grad = objective_vg(jnp.asarray(log_w_np))
        return float(val), np.asarray(grad, dtype=np.float64)

    opt_res = minimize(
        fun=scipy_objective,
        x0=log_width_init,
        method="BFGS",
        jac=True,
        options={"maxiter": maxiter},
    )

    # Final evaluation
    var_opt, e_opt = compute_variance_and_energy(opt_res.x, n)
    widths_opt = np.exp(np.asarray(opt_res.x))

    return {
        "n": n,
        "variance": var_opt,
        "energy": e_opt,
        "widths_opt": widths_opt,
        "success": bool(opt_res.success),
        "message": opt_res.message,
        "niter": opt_res.nit,
    }


def optimize_logspace_for_n(
    n: int,
    bounds_init: np.ndarray | None = None,
    objective: str = "variance",
    maxiter: int = 500,
) -> dict:
    """Optimize only minval and maxval for log-spaced Gaussian widths.

    The n widths are determined by np.log(np.logspace(minval, maxval, n)),
    so widths range from 10^minval to 10^maxval in log-spacing.

    Args:
        n: Basis size (number of Gaussians)
        bounds_init: Initial [minval, maxval]. If None, use [-1, 2] (widths from 0.1 to 100).
        objective: What to minimize - "energy" or "variance"
        maxiter: Maximum optimization iterations

    Returns:
        Dictionary with optimization results including variance, energy, and optimal widths.
    """
    if bounds_init is None:
        bounds_init = np.array([-1.0, 2.0])

    # Select objective function
    if objective == "energy":
        objective_vg = make_logspace_energy_objective(n)
    elif objective == "variance":
        objective_vg = make_logspace_variance_objective(n)
    else:
        raise ValueError(f"Unknown objective: {objective}. Use 'energy' or 'variance'.")

    def scipy_objective(bounds_np):
        val, grad = objective_vg(jnp.asarray(bounds_np))
        return float(val), np.asarray(grad, dtype=np.float64)

    opt_res = minimize(
        fun=scipy_objective,
        x0=bounds_init,
        method="BFGS",
        jac=True,
        options={"maxiter": maxiter, "gtol": 1e-16},
    )

    # Final evaluation
    minval_opt, maxval_opt = opt_res.x
    log_w_opt = np.log(np.logspace(minval_opt, maxval_opt, n))
    var_opt, e_opt = compute_variance_and_energy(log_w_opt, n)
    widths_opt = np.exp(log_w_opt)

    return {
        "n": n,
        "variance": var_opt,
        "energy": e_opt,
        "widths_opt": widths_opt,
        "minval": minval_opt,
        "maxval": maxval_opt,
        "success": bool(opt_res.success),
        "message": opt_res.message,
        "niter": opt_res.nit,
    }


# =============================================================================
# Helper for building initial guesses from previous solutions
# =============================================================================


def expand_log_widths(prev_logs: np.ndarray, spread_factor: float = 0.05) -> np.ndarray:
    """Expand a set of log-widths by inserting one in the middle and spreading the rest.

    Example: [2, 5, 7, 10] with 5% spread -> [1.9, 4.92, 6, 7.12, 10.5]
    """
    prev_sorted = np.sort(prev_logs)
    n_prev = len(prev_sorted)

    # Insert new element at midpoint between two middle elements
    insert_idx = n_prev // 2
    middle_log = 0.5 * (prev_sorted[insert_idx - 1] + prev_sorted[insert_idx])

    # Pull apart: linear shift from -spread_factor to +spread_factor
    shifts = np.linspace(-spread_factor, spread_factor, n_prev)
    spread_logs = prev_sorted * (1 + shifts)

    return np.insert(spread_logs, insert_idx, middle_log)


# =============================================================================
# Main entry point
# =============================================================================


def main():
    # Known good starting solution for n=6
    n_6_widths = np.array(
        [
            0.42165816943491824,
            0.8728668410863247,
            1.8629775886223578,
            4.211462150107637,
            10.316307838808878,
            28.30265815953873,
        ]
    )

    results = []
    prev_logs = np.log(n_6_widths)
    bounds_init = None
    for n in range(15, 30):
        log_init = expand_log_widths(prev_logs, spread_factor=0.05)
        res = optimize_logspace_for_n(n, bounds_init=bounds_init, objective="variance")
        results.append(res)

        print(
            f"n={res['n']:2d} | success={res['success']} | "
            f"var={res['variance']:.3e} | E={res['energy']:.6f} | iters={res['niter']}"
        )
        print(f"n={res['n']:2d}: widths={list(res['widths_opt'])}")

        prev_logs = np.log(res["widths_opt"])
        bounds_init = np.log10(np.array([min(res["widths_opt"]), max(res["widths_opt"])]))
        print("Bounds init for next n:", bounds_init)


if __name__ == "__main__":
    main()
