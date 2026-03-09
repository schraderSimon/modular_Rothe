"""Propagate a single Gaussian in a 1D symmetric double-well potential.

Potential: V(x) = x^4/4 - x^2/2
  - Two minima at x = ±1, V(±1) = -1/4
  - Barrier at x = 0, height 1/4 above the minima

The initial state is a single Gaussian centred at the right minimum (x=+1)
with zero momentum and width matched to the local harmonic approximation
(V''(±1) = 2, so ω = √2, giving a = ω^{1/2}/√2 = 2^{1/4}/2 ≈ 0.59).

This is the simplest possible sanity-check: with a single Gaussian the code
cannot represent tunnelling faithfully (you'd need ≥ 2 Gaussians), but you
can watch the code add basis functions as the Rothe error budget forces it.

Usage:
    # minimal run — 200 steps of dt=0.1 (t_final=20), generous epsilon
    python propagate_double_well.py --epsilon 1e-3

    # tighter tolerance will trigger Gaussian addition
    python propagate_double_well.py --epsilon 1e-5 --dt 0.05

    # start with more than one Gaussian
    python propagate_double_well.py --epsilon 1e-4 --num_gaussians 4
"""

import argparse

import jax
import jax.numpy as jnp

jax.config.update("jax_compilation_cache_dir", str(__import__("pathlib").Path.home() / ".cache/jax"))
jax.config.update("jax_persistent_cache_min_compile_time_secs", 1.0)

from rothe.io import OutputConfig
from rothe.solver import RotheSolver, setUpRotheErrorAndGradient_jit
from rothe.systems.common import eval_initial_state_properties, resume_or_initialize_rothe, set_up_SHH2

# ---------------------------------------------------------------------------
# Potential
# ---------------------------------------------------------------------------

DOUBLE_WELL_STRING = """
dimension 1
polynomial
x0x0x0x0: 0.25
x0x0: -0.5
"""

# ---------------------------------------------------------------------------


def make_initial_params(n: int, width: float = 0.59) -> jnp.ndarray:
    """Place `n` Gaussians near the right minimum (x=+1) with slightly varied centres.

    With n=1 the single Gaussian sits exactly at the minimum.  With n>1 the
    additional ones are spread symmetrically around it so the optimiser has
    something to work with from the first step.

    Parameters
    ----------
    n : int       Number of Gaussians.
    width : float Gaussian width parameter `a` (default matches local ω = √2).
    """
    params = jnp.zeros((n, 1, 4), dtype=jnp.float64)
    params = params.at[:, :, 0].set(width)  # a  (width)
    # b=0, p=0 are already zero — no initial chirp / momentum

    # centres: spread evenly from 0.5 to 1.5 so they don't all pile on top
    centres = jnp.linspace(0.5, 1.5, n) if n > 1 else jnp.array([1.0])
    params = params.at[:, 0, 2].set(centres)  # mu
    return params


def parse_args():
    parser = argparse.ArgumentParser(description="Propagate a single Gaussian in a 1D double-well.")
    parser.add_argument("--num_gaussians", type=int, default=1, help="Initial number of Gaussians")
    parser.add_argument("--dt", type=float, default=0.1, help="Time step (default 0.1)")
    parser.add_argument("--t_final", type=float, default=20.0, help="Final time (default 20)")
    parser.add_argument("--epsilon", type=float, default=1e-3, help="Global Rothe-error budget (default 1e-3)")
    parser.add_argument("--splitting_type", type=str, default="none", choices=["none", "kinetic"])
    parser.add_argument("--random_seed", type=int, default=0)
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--t_start", type=float, default=None, help="Resume from checkpoint at or before this time")
    return parser.parse_args()


def main():
    args = parse_args()
    dt = args.dt
    epsilon = args.epsilon
    n = args.num_gaussians
    splitting_type = args.splitting_type
    random_seed = args.random_seed
    nsteps_total = int(round(args.t_final / dt))

    SHH2 = set_up_SHH2(potential_string=DOUBLE_WELL_STRING, D=1)
    rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type)

    params_init = make_initial_params(n)
    coeffs_init = jnp.ones(n, dtype=jnp.complex128)
    S0, _, _ = SHH2(0.0, params_init, params_init, splitting_type=splitting_type)
    coeffs_init = coeffs_init / jnp.sqrt(jnp.vdot(coeffs_init, S0 @ coeffs_init))

    print("=== Initial state ===")
    eval_initial_state_properties(coeffs_init, params_init, SHH2)

    initial_error = rothe_error(params_init, params_init, coeffs_init, SHH2, 0.0, dt)
    print(f"Initial Rothe error: {float(initial_error):.3e}")

    sim_name = f"DoubleWell_1D_n{n}_eps{epsilon:.0e}"

    rothesolver = RotheSolver(
        SHH2=SHH2,
        dt=dt,
        t=0.0,
        epsilon=epsilon,
        random_seed=random_seed,
        params_old=params_init,
        coeffs_old=coeffs_init,
        rothe_grad_fn=rothe_vg_jit,
        rothe_nograd=rothe_error,
        splitting_type=splitting_type,
        output_config=OutputConfig(name=sim_name, polynomial_string=DOUBLE_WELL_STRING, epsilon=epsilon),
    )

    nsteps = resume_or_initialize_rothe(
        rothesolver=rothesolver,
        nsteps_total=nsteps_total,
        solver_name=sim_name,
        splitting_type=splitting_type,
        potential_string=DOUBLE_WELL_STRING,
        epsilon=epsilon,
        dt=dt,
        initial_error=initial_error,
        params=params_init,
        coeffs=coeffs_init,
        t_start=args.t_start,
    )

    rothesolver.propagate(nsteps, num_time_steps_total=nsteps_total, maxiter=args.maxiter)


if __name__ == "__main__":
    main()
