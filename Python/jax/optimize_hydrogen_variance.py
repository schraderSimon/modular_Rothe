import numpy as np
import pandas as pd
from scipy.optimize import minimize

import jax
import jax.numpy as jnp
from general_wf import generalPotentialSolver
from utils import calculate_variance

# Build the hydrogen potential string from stored Gaussian fit coefficients
# The linear coefficients are negated to represent the attractive Coulomb tail.
df = pd.read_csv("gaussian_Coulomb/coeffs_mu10_N20.csv")
lincoeffs = df["linear"]
widths = df["nonlinear"]

hydrogen_string = """
    dimension 3
    exponential
    """
for lc, w in zip(lincoeffs, widths):
    hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"


def optimize_variance_for_n(n: int, log_width_init: np.ndarray | None = None):
    """Minimize ground-state variance over isotropic Gaussian widths for a basis of size n.

    Args:
        n: basis size
        log_width_init: optional initial guess (log widths). If None, use a spread grid.
    """

    # Initial log-widths span ~[0.1, 100] if no prior guess is given
    if log_width_init is None:
        log_width_init = np.logspace(-1.0, 1, n)

    def params_from_log_widths(log_w):
        w = jnp.exp(jnp.asarray(log_w))  # ensure positivity
        params = jnp.zeros((n, 3, 4), dtype=jnp.float64)
        return params.at[:, :, 0].set(w[:, None])

    @jax.jit
    def variance_only(log_w):
        # Build a fresh solver per call to keep the function side-effect free for autodiff.
        params = params_from_log_widths(log_w)
        system = generalPotentialSolver(params, hydrogen_string)
        S = system.calculate_S()
        H = system.calculate_H()
        H2 = system.calculate_H2()

        # JAX-friendly orthogonalization without boolean masking (avoid NonConcreteBooleanIndexError).
        w, U = jnp.linalg.eigh(0.5 * (S + S.conj().T))
        w_clipped = jnp.clip(jnp.real(w), 1e-12, None)
        X = U * (1.0 / jnp.sqrt(w_clipped))[None, :]
        H_orth = X.conj().T @ H @ X
        H2_orth = X.conj().T @ H2 @ X

        eigvals, eigvecs_orth = jnp.linalg.eigh(H_orth)
        energy = jnp.real(eigvals[0])
        gs = eigvecs_orth[:, 0]
        variance = calculate_variance(gs, H2_orth, energy)
        return variance

    def energy_only(log_w):
        params = params_from_log_widths(log_w)
        system = generalPotentialSolver(params, hydrogen_string)
        S = system.calculate_S()
        H = system.calculate_H()
        H2 = system.calculate_H2()

        # JAX-friendly orthogonalization without boolean masking (avoid NonConcreteBooleanIndexError).
        w, U = jnp.linalg.eigh(0.5 * (S + S.conj().T))
        w_clipped = jnp.clip(jnp.real(w), 1e-16, None)
        X = U * (1.0 / jnp.sqrt(w_clipped))[None, :]
        H_orth = X.conj().T @ H @ X
        H2_orth = X.conj().T @ H2 @ X

        eigvals, eigvecs_orth = jnp.linalg.eigh(H_orth)
        energy = jnp.real(eigvals[0])
        return energy

    # Jitted value and grad of the variance
    # variance_vg = jax.jit(jax.value_and_grad(variance_only))
    energy_vg = jax.jit(jax.value_and_grad(energy_only))

    def objective_with_grad(log_w_np):
        val, grad = energy_vg(jnp.asarray(log_w_np))
        return float(val), np.asarray(grad, dtype=np.float64)

    opt_res = minimize(
        fun=objective_with_grad,
        x0=log_width_init,
        method="BFGS",
        jac=True,
        options={"maxiter": 500},
    )

    # Final evaluation (variance and energy) using the optimized widths
    def variance_and_energy(log_w):
        params = params_from_log_widths(log_w)
        system = generalPotentialSolver(params, hydrogen_string)
        S = system.calculate_S()
        H = system.calculate_H()
        H2 = system.calculate_H2()

        w, U = jnp.linalg.eigh(0.5 * (S + S.conj().T))
        w_clipped = jnp.clip(jnp.real(w), 1e-12, None)
        X = U * (1.0 / jnp.sqrt(w_clipped))[None, :]
        H_orth = X.conj().T @ H @ X
        H2_orth = X.conj().T @ H2 @ X
        eigvals, eigvecs_orth = jnp.linalg.eigh(H_orth)
        energy = float(jnp.real(eigvals[0]))
        gs = eigvecs_orth[:, 0]
        variance = float(calculate_variance(gs, H2_orth, energy))
        return variance, energy

    var_opt, e_opt = variance_and_energy(opt_res.x)
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


def main():
    results = []
    prev_logs = None
    for n in range(2, 21):
        if prev_logs is None:
            log_init = None
        else:
            # Insert one extra Gaussian roughly in the middle of the span of previous widths
            prev_sorted = np.sort(prev_logs)
            middle_log = np.mean(prev_sorted)
            insert_idx = len(prev_sorted) // 2
            log_init = np.insert(prev_sorted, insert_idx, middle_log)
            log_init = log_init * np.random.uniform(
                0.99, 1.01, size=log_init.shape
            )  # slight random perturbation
        log_init = None
        res = optimize_variance_for_n(n, log_width_init=log_init)
        results.append(res)
        print(
            f"n={res['n']:2d} | success={res['success']} | var={res['variance']:.3e} | E={res['energy']:.6f} | iters={res['niter']}"
        )
        res = results[-1]
        log_w = np.log(res["widths_opt"])
        print(f"n={res['n']:2d}: widths={list(np.array(res['widths_opt']))}")
        prev_logs = log_w


if __name__ == "__main__":
    main()
