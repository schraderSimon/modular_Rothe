"""Propagate Gaussians in a 1D symmetric double-well: V(x) = x^4/4 - x^2/2."""

import argparse

import jax
import jax.numpy as jnp

jax.config.update("jax_compilation_cache_dir", str(__import__("pathlib").Path.home() / ".cache/jax"))
jax.config.update("jax_persistent_cache_min_compile_time_secs", 1.0)

from rothe.io import OutputConfig
from rothe.solver import RotheSolver, setUpRotheErrorAndGradient_jit
from rothe.systems.common import eval_initial_state_properties, resume_or_initialize_rothe, set_up_SHH2

DOUBLE_WELL_STRING = """
dimension 1
polynomial
x0x0x0x0: 0.25
x0x0: -0.5
"""

RIGHT_WELL_CENTER = 1
LEFT_WELL_CENTER = -1
GROUND_STATE_WIDTH = 1 / 2 ** (1 / 4)


def make_initial_params(n: int, width: float = GROUND_STATE_WIDTH) -> jnp.ndarray:
    params = jnp.zeros((n, 1, 4), dtype=jnp.float64)
    params = params.at[:, :, 0].set(width)
    if n == 1:
        centers = jnp.array([RIGHT_WELL_CENTER])
    else:
        centers = jnp.linspace(LEFT_WELL_CENTER, RIGHT_WELL_CENTER, n, dtype=jnp.float64)
    params = params.at[:, 0, 2].set(centers)
    return params


def make_initial_coeffs(n: int) -> jnp.ndarray:
    coeffs = jnp.zeros(n, dtype=jnp.complex128)
    return coeffs.at[0].set(1.0 + 0.0j)


def parse_args():
    parser = argparse.ArgumentParser(description="Propagate a single Gaussian in a 1D double-well.")
    parser.add_argument("--num_gaussians", type=int, default=25, help="Initial number of Gaussians")
    parser.add_argument("--dt", type=float, default=0.1, help="Time step (default 0.1)")
    parser.add_argument("--t_final", type=float, default=20.0, help="Final time (default 20)")
    parser.add_argument("--epsilon", type=float, default=1e-2, help="Global Rothe-error budget (default 1e-3)")
    parser.add_argument("--regularization_lambda", type=float, default=1e-10, help="Tikhonov regularization strength")
    parser.add_argument("--splitting_type", type=str, default="none", choices=["none", "kinetic"])
    parser.add_argument("--random_seed", type=int, default=0)
    parser.add_argument("--maxiter", type=int, default=200)
    parser.add_argument("--t_start", type=float, default=None, help="Resume from checkpoint at or before this time")
    parser.add_argument(
        "--remove_overlapping_gaussians",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Remove one dynamic Gaussian if its overlap with a frozen Gaussian exceeds 0.99, then re-optimize.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    dt = args.dt
    epsilon = args.epsilon
    regularization_lambda = args.regularization_lambda
    n = args.num_gaussians
    splitting_type = args.splitting_type
    random_seed = args.random_seed
    nsteps_total = int(round(args.t_final / dt))

    SHH2 = set_up_SHH2(potential_string=DOUBLE_WELL_STRING, D=1)
    rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type)

    params_init = make_initial_params(n)
    coeffs_init = make_initial_coeffs(n)
    eval_initial_state_properties(coeffs_init, params_init, SHH2)

    initial_error = rothe_error(params_init, params_init, coeffs_init, SHH2, 0.0, dt, lambda_=regularization_lambda)
    print(f"Initial Rothe error: {float(initial_error):.3e}")

    sim_name = f"DoubleWell_1D_n{n}_eps{epsilon:.0e}_lambda{regularization_lambda:.0e}"

    rothesolver = RotheSolver(
        SHH2=SHH2,
        dt=dt,
        t=0.0,
        epsilon=epsilon,
        regularization_lambda=regularization_lambda,
        random_seed=random_seed,
        params_old=params_init,
        coeffs_old=coeffs_init,
        rothe_grad_fn=rothe_vg_jit,
        rothe_nograd=rothe_error,
        splitting_type=splitting_type,
        output_config=OutputConfig(name=sim_name, polynomial_string=DOUBLE_WELL_STRING, epsilon=epsilon),
        remove_overlapping_gaussians=args.remove_overlapping_gaussians,
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
