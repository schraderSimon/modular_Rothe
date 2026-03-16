"""Propagate the hydrogen atom wave function using the Rothe method."""

import argparse

import numpy as np

import jax
import jax.numpy as jnp

jax.config.update("jax_compilation_cache_dir", str(__import__("pathlib").Path.home() / ".cache/jax"))
jax.config.update("jax_persistent_cache_min_compile_time_secs", 1.0)
from rothe.io import OutputConfig
from rothe.solver import RotheSolver, setUpRotheErrorAndGradient_jit
from rothe.systems.common import ExternalFieldParams, resume_or_initialize_rothe, set_up_SHH2
from rothe.systems.hydrogen import (
    make_hydrogen_squared_string,
    make_hydrogen_string,
    read_best_coefficients,
    set_up_hydrogen_string_name,
)


def make_dynamic_params(n_per_side=10, width=1 / np.sqrt(2), D=4, z_spacing=0.5):
    """Create 2*n_per_side Gaussians placed symmetrically along the z-axis."""
    n_total = 2 * n_per_side
    params = jnp.zeros((n_total, D, 4), dtype=jnp.float64)
    params = params.at[:, :, 0].set(width)

    z_positions_pos = jnp.arange(1, n_per_side + 1) * z_spacing
    z_positions_neg = -z_positions_pos
    z_positions = jnp.concatenate([z_positions_pos, z_positions_neg])  # (2*n_per_side,)
    params = params.at[:, 2, 2].set(z_positions)  # dim 2 (z), column 2 (mu)
    return params


def parse_args():
    parser = argparse.ArgumentParser(description="Propagate the hydrogen atom wave function using the Rothe method.")
    parser.add_argument("--num_gauss_potential", type=int, default=21)
    parser.add_argument("--num_gauss_wavefunction", type=int)
    parser.add_argument("--dt", type=complex, help="time step", default=0.2)
    parser.add_argument(
        "--splitting_type",
        type=str,
        default="none",
        choices=["none", "kinetic"],
        help="Rothe splitting mode to benchmark",
    )
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--regularization_lambda", type=float, default=1e-10, help="Tikhonov regularization strength")
    parser.add_argument("--random_seed", type=int, default=0, help="RNG seed for Gaussian sampling")
    parser.add_argument("--maxiter", type=int, default=50)
    parser.add_argument("--num_dynamic", type=int, default=10, help="Number of dynamic Gaussians per side")
    parser.add_argument(
        "--t_start",
        type=float,
        default=None,
        help="Time to resume propagation from (resumes from checkpoint at or before this time)",
    )
    parser.add_argument(
        "--num_gauss_fit",
        type=int,
        default=None,
        help="Number of Gaussians in the pre-fitted V^2 approximation. "
        "If given, uses square_coeffs CSV instead of exact squaring.",
    )
    parser.add_argument(
        "--frozen",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Freeze ground-state Gaussians (default: True). " "Use --no-frozen to make all Gaussians optimizable.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    splitting_type = args.splitting_type
    ngp = num_gauss_potential = args.num_gauss_potential
    ngwf = num_gauss_wavefunction = args.num_gauss_wavefunction
    maxiter = args.maxiter
    dt = args.dt
    epsilon = args.epsilon
    regularization_lambda = args.regularization_lambda
    random_seed = args.random_seed
    t_start = args.t_start
    is_imaginary_time = bool(dt.imag)
    if not is_imaginary_time:
        dt = float(dt.real)

    external_field_params = ExternalFieldParams()
    if not is_imaginary_time:
        external_field_params.E0 = 0.06
    else:
        external_field_params.E0 = 0.0
    use_frozen = args.frozen and not is_imaginary_time
    num_dynamic_eff = 0 if (is_imaginary_time or not use_frozen) else args.num_dynamic
    string_name = set_up_hydrogen_string_name(
        external_field_params, ngp, ngwf, dt, epsilon, regularization_lambda, num_dynamic=num_dynamic_eff
    )
    tfinal = external_field_params.t_final
    nsteps = int(jnp.ceil(tfinal / jnp.abs(dt)))

    # Set up the initial wavefunction and potential
    hydrogen_string = make_hydrogen_string(num_gauss_potential, external_field_params)
    try:
        if dt.imag:
            raise FileNotFoundError(
                "Imaginary time propagation should start from random initial state, not ground state data."
            )
        initial_linear_coeffs, params_init = read_best_coefficients(
            num_gauss_wavefunction=num_gauss_wavefunction, num_gauss_potential=num_gauss_potential
        )
    except FileNotFoundError:
        if dt.imag:
            print(
                f"No ground state data for ng_wf={ngwf}, ng_pot={ngp}. "
                f"Using random initial state for imaginary time propagation."
            )
            D = 3  # Hydrogen is 3D
            widths = np.logspace(-2, 2, num_gauss_wavefunction)
            params_init = jnp.zeros((num_gauss_wavefunction, D, 4), dtype=jnp.float64)
            params_init = params_init.at[:, :, 0].set(widths[:, None])
            initial_linear_coeffs = jnp.ones(num_gauss_wavefunction, dtype=jnp.complex128)
        else:
            raise FileNotFoundError(
                f"No ground state data for ng_wf={ngwf}, ng_pot={ngp}. "
                f"Cannot proceed with real time propagation without initial state."
            )
    squared_string = (
        make_hydrogen_squared_string(num_gauss_potential, args.num_gauss_fit)
        if args.num_gauss_fit is not None
        else None
    )

    if not use_frozen:
        # All Gaussians are dynamic, nothing is frozen.
        params_all = jnp.asarray(params_init)
        coeffs_all = jnp.asarray(initial_linear_coeffs)
        params_frozen_solver = None
    else:
        # Real-time propagation with frozen core: ground-state Gaussians are
        # frozen; prepend dynamic Gaussians that are free to move.
        params_frozen_solver = params_init
        params_dynamic = make_dynamic_params(D=int(params_init.shape[1]), n_per_side=args.num_dynamic)
        n_dynamic = params_dynamic.shape[0]
        params_all = jnp.concatenate([params_dynamic, params_frozen_solver], axis=0)
        # Dynamic part starts with zero coefficients; frozen part carries the
        # ground-state coefficients.
        coeffs_all = jnp.concatenate([jnp.zeros(n_dynamic, dtype=jnp.complex128), jnp.asarray(initial_linear_coeffs)])

    # Set up SHH2 and Rothe error/gradient functions
    SHH2 = set_up_SHH2(potential_string=hydrogen_string, D=3, squared_string=squared_string)
    rothe_error, rothe_gradient = setUpRotheErrorAndGradient_jit(splitting_type)

    # Set up Rothe solver
    rothesolver = RotheSolver(
        SHH2=SHH2,
        dt=dt,
        t=0,
        epsilon=epsilon,
        regularization_lambda=regularization_lambda,
        random_seed=random_seed,
        params_old=params_all,
        coeffs_old=coeffs_all,
        rothe_grad_fn=rothe_gradient,
        rothe_nograd=rothe_error,
        splitting_type=splitting_type,
        output_config=OutputConfig(name=string_name, polynomial_string=hydrogen_string, epsilon=epsilon),
        params_frozen=params_frozen_solver,
    )
    nsteps_total = nsteps
    nsteps = resume_or_initialize_rothe(
        rothesolver=rothesolver,
        nsteps_total=nsteps,
        solver_name=rothesolver.output_config.name,
        splitting_type=splitting_type,
        potential_string=hydrogen_string,
        epsilon=epsilon,
        dt=dt,
        initial_error=0,
        params=params_all,
        coeffs=coeffs_all,
        t_start=t_start,
    )

    rothesolver.propagate(nsteps, num_time_steps_total=nsteps_total, maxiter=maxiter)


if __name__ == "__main__":
    main()
