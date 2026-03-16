"""Simulate the Henon-Heiles system using the Rothe method."""

import argparse

import jax
import jax.numpy as jnp
from rothe.io import OutputConfig
from rothe.solver import RotheSolver, setUpRotheErrorAndGradient_jit
from rothe.systems.common import resume_or_initialize_rothe, set_up_SHH2
from rothe.systems.henon_heiles import initialize_henon_heiles_params, make_henon_heiles_string

jax.config.update("jax_compilation_cache_dir", str(__import__("pathlib").Path.home() / ".cache/jax"))


def parse_args():
    """Parse command-line arguments for Henon-Heiles simulation."""
    parser = argparse.ArgumentParser(description="Simulate Henon-Heiles system using the Rothe method.")
    parser.add_argument("--dim", type=int, default=2, help="Dimension of the system (1, 2, 4, or 6)")
    parser.add_argument("--splitting_type", type=str, default="none", help="Splitting type for Rothe propagation")
    parser.add_argument("--dt", type=float, default=0.01, help="Time step")
    parser.add_argument("--num_gaussians", type=int, default=1, help="Number of Gaussians")
    parser.add_argument("--epsilon", type=float, required=True, help="Global Rothe-error budget")
    parser.add_argument("--regularization_lambda", type=float, default=1e-10, help="Tikhonov regularization strength")
    parser.add_argument("--random_seed", type=int, default=0, help="RNG seed for Gaussian sampling")
    return parser.parse_args()


def main():
    args = parse_args()
    splitting_type = args.splitting_type
    D = args.dim
    dt = args.dt
    n = args.num_gaussians
    epsilon = args.epsilon
    regularization_lambda = args.regularization_lambda
    random_seed = args.random_seed
    nsteps = int(100 / dt)

    # System configuration
    hehe_string = make_henon_heiles_string(D)

    # Initialize parameters
    params_init = initialize_henon_heiles_params(n=n, D=D, seed=42)

    # Set up SHH2 and Rothe error/gradient functions
    SHH2 = set_up_SHH2(potential_string=hehe_string, D=D)
    rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type)

    # Set up wavefunction coefficients
    coeffs_init = jnp.zeros(n, dtype=jnp.complex128)
    coeffs_init = coeffs_init.at[0].set(1.0)
    S0, _, _ = SHH2(0, params_init, params_init, splitting_type=splitting_type)
    coeffs_init = coeffs_init / jnp.sqrt(jnp.vdot(coeffs_init, S0 @ coeffs_init))

    filename = f"HenonHeiles_{D}D_{n}__epsilon={epsilon:.1e}__lambda={regularization_lambda:.1e}"
    # Set up Rothe solver
    rothesolver = RotheSolver(
        SHH2=SHH2,
        dt=dt,
        t=0,
        epsilon=epsilon,
        regularization_lambda=regularization_lambda,
        random_seed=random_seed,
        params_old=params_init,
        coeffs_old=coeffs_init,
        rothe_grad_fn=rothe_vg_jit,
        rothe_nograd=rothe_error,
        splitting_type=splitting_type,
        output_config=OutputConfig(name=filename, polynomial_string=hehe_string, epsilon=epsilon),
    )

    nsteps_total = nsteps
    nsteps = resume_or_initialize_rothe(
        rothesolver=rothesolver,
        nsteps_total=nsteps,
        solver_name=rothesolver.output_config.name,
        splitting_type=splitting_type,
        potential_string=hehe_string,
        epsilon=epsilon,
        dt=dt,
        initial_error=0,
        params=params_init,
        coeffs=coeffs_init,
    )

    rothesolver.propagate(nsteps, num_time_steps_total=nsteps_total, maxiter=1000)


if __name__ == "__main__":
    main()
