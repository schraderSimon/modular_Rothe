import os
import sys

import numpy as np

import jax.numpy as jnp
from general_wf import generalPotentialSolver
from jax_Rothe import RotheSolver, save_rothe_state, setUpRotheErrorAndGradient_jit
from read_string import read_string

splitting_type = sys.argv[1] if len(sys.argv) > 1 else "none"
# Henon-Heiles potential
example_string = """
    dimension 6
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x2x2: 0.5
    x3x3: 0.5
    x4x4: 0.5
    x5x5: 0.5
    x0x0x1: 0.111803
    x1x1x2: 0.111803
    x2x3x3: 0.111803
    x3x4x4: 0.111803
    x4x5x5: 0.111803
    x1x1x1: -0.03726766666
    x2x2x2: -0.03726766666
    x3x3x3: -0.03726766666
    x4x4x4: -0.03726766666
    x5x5x5: -0.03726766666
    """
example_string = """
    dimension 4
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x2x2: 0.5
    x3x3: 0.5
    x0x0x1: 0.111803
    x1x1x2: 0.111803
    x2x3x3: 0.111803
    x1x1x1: -0.03726766666
    x2x2x2: -0.03726766666
    x3x3x3: -0.03726766666
    """

polynomial = read_string(example_string)
rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type=splitting_type)
dt = 0.01
n = 10
key, value = next(iter(polynomial.items()))
D = len(key)


def make_SHH2_pure(splitting_type="none"):
    osc = generalPotentialSolver(p_init, example_string)

    def SHH2(t, params, splitting_type=splitting_type):
        osc.update_parameters(params)
        return osc.calculate_SHH2(t=t, splitting_type=splitting_type)

    return SHH2


p_init = jnp.zeros((n, D, 4), dtype=jnp.float64)
p_init = p_init.at[:, :, 0].set(1 / jnp.sqrt(2))  # Width parameters
p_init = p_init.at[0, :, 2].set(2.0)  # Set mu to (2,...,2)
for i in range(1, n):
    p_init = p_init.at[i, :, 2].set(
        2 + np.random.uniform(-0.5, 0.5, (D,))
    )  # Move to the right
# ---------- build WF + solver ---------------------------------------
oscillator = generalPotentialSolver(p_init, example_string)
SHH2 = make_SHH2_pure(splitting_type=splitting_type)  # The function to calculate S, H, H2
Smat = oscillator.calculate_S()
c_init = jnp.array([1.0] * n, dtype=jnp.complex128)  # c
c_init = jnp.zeros_like(c_init)
c_init = c_init.at[0].set(1.0)
norm = c_init.T @ Smat @ c_init
c_init = c_init / jnp.sqrt(norm)
initial_error = rothe_error(p_init, p_init, c_init, SHH2, 0, dt)
rothesolver = RotheSolver(
    SHH2,
    dt,
    0,
    p_init,
    c_init,
    rothe_grad_jax=rothe_vg_jit,
    rothe_nograd=rothe_error,
    splitting_type=splitting_type,
    name="HenonHeiles_%dD_%d" % (D, n),
    polynomial_string=example_string,
)
nsteps = int(100 / dt)
try:
    info = rothesolver.resume_from_file(t_start=None)  # or pick a time
    print(f"Resumed at t={info['t']} (step {info['idx']}), trimmed={info['trimmed']}")
    nsteps = nsteps - info["idx"]  # remaining steps
except FileNotFoundError:
    # optional: log initial state at t=0, so resumes have an anchor
    save_rothe_state(
        rothesolver.name,
        splitting_type,
        example_string,
        dt,
        0.0,
        float(initial_error),
        p_init,
        c_init,
    )
    print("No prior file. Starting fresh at t=0.")

rothesolver.propagate(nsteps, maxiter=1000)
