"""
Common system helpers shared across physical systems.
"""

from attr import dataclass

import jax.numpy as jnp
from rothe.wavefunction import generalPotentialSolver


@dataclass
class ExternalFieldParams:
    """Parameters for a time-dependent laser field plus envelope."""

    E0: float = 0.06  # field strength in atomic units
    omega: float = 0.057  # frequency in atomic units
    N_T: int = 3  # Number of cycles

    @property
    def intensity(self) -> float:
        return 3.5094e16 * self.E0**2

    @property
    def t_final(self) -> float:
        return 2 * jnp.pi / self.omega * self.N_T


def set_up_SHH2(potential_string, D, squared_string=None):
    """Build the SHH2 function for a given potential.
    Input:
    - potential_string: string specifying the potential
    - D: spatial dimension
    - squared_string: string specifying the squared potential (optional)
    Returns:
        SHH2: function(t, params_ket, params_bra, splitting_type) -> (S, H, H2)
    """
    params_dummy = jnp.ones((1, D, 4))
    system = generalPotentialSolver(params_dummy, params_dummy, potential_string, squared_string=squared_string)

    def SHH2(t, params_ket, params_bra, splitting_type="none"):
        system.update_parameters(params_ket, params_bra)
        return system.calculate_SHH2(t=t, splitting_type=splitting_type)

    return SHH2


def eval_initial_state_properties(initial_linear_coeffs, params_init, SHH2):
    """Evaluate and print norm, energy, and variance of an initial state."""
    initial_linear_coeffs = jnp.array(initial_linear_coeffs)
    S, H, H2 = SHH2(0, params_init, params_init, "none")
    norm = jnp.vdot(initial_linear_coeffs, S @ initial_linear_coeffs)
    initial_linear_coeffs = initial_linear_coeffs / jnp.sqrt(norm)
    E0 = jnp.vdot(initial_linear_coeffs, H @ initial_linear_coeffs)
    Esq = jnp.vdot(initial_linear_coeffs, H2 @ initial_linear_coeffs)
    print("Initial energy:", E0)
    print("Initial norm:", norm)
    print("Initial variance:", Esq - E0**2)


def resume_or_initialize_rothe(
    rothesolver,
    nsteps_total,
    solver_name,
    splitting_type,
    potential_string,
    epsilon,
    dt,
    initial_error=0.0,
    params=None,
    coeffs=None,
    t_start=None,
):
    """
    Resume Rothe solver from checkpoint or initialize from scratch.

    Returns:
        nsteps_remaining: number of steps left to propagate
    """
    from ..io import save_rothe_state

    try:
        info = rothesolver.resume_from_file(t_start=t_start)
        print(f"Resumed at t={info['t']} (step {info['idx']}), trimmed={info['trimmed']}")
        return nsteps_total - info["idx"]
    except FileNotFoundError:
        n_frozen = rothesolver.params_frozen.shape[0] if rothesolver.params_frozen is not None else 0
        n_dynamic_init = params.shape[0] - n_frozen if params is not None else 0
        rng_state = rothesolver._rng.bit_generator.state if hasattr(rothesolver, "_rng") else None
        save_rothe_state(
            solver_name,
            splitting_type,
            potential_string,
            epsilon,
            dt,
            0.0,
            float(initial_error),
            params=params,
            coeffs=coeffs,
            n_dynamic=n_dynamic_init,
            rng_state=rng_state,
            path=rothesolver.out_dir,
        )
        print("No prior file. Starting fresh at t=0.")
        return nsteps_total
