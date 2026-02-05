import os
import sys

os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true --xla_gpu_autotune_level=4"

import numpy as np

import jax.numpy as jnp
from libraries.general_wf import generalPotentialSolver
from libraries.jax_Rothe import RotheSolver, save_rothe_state, setUpRotheErrorAndGradient_jit
from libraries.read_string import read_string

splitting_type = sys.argv[1] if len(sys.argv) > 1 else "none"
# Henon-Heiles potential
example_string_6D = """
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

example_string_4D = """
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
example_string_2D = """
    dimension 2
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x0x0x1: 0.111803
    x1x1x1: -0.03726766666
    x0x0x0x0: 0.0007812444255625
    x1x1x1x1: 0.0007812444255625
    x0x0x1x1: 0.001562488851125
    """
example_string_1D = """
    dimension 1
    polynomial
    x0x0: 0.5
    x0x0x0x0: 0.0007812444255625
    """
dim = 2
if dim == 6:
    example_string = example_string_6D
elif dim == 4:
    example_string = example_string_4D
elif dim == 1:
    example_string = example_string_1D
else:
    example_string = example_string_2D
polynomial, exponential = read_string(example_string)
rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type=splitting_type)
dt = 0.01
n = 1
key, value = next(iter(polynomial.items()))
D = len(key)


def make_SHH2_pure(splitting_type="none"):
    osc = generalPotentialSolver(params_init, example_string)

    def SHH2(t, params, splitting_type=splitting_type):
        osc.update_parameters(params)
        return osc.calculate_SHH2(t=t, splitting_type=splitting_type)

    return SHH2


np.random.seed(42)

params_init = jnp.zeros((n, D, 4), dtype=jnp.float64)
params_init = params_init.at[:, :, 0].set(1 / jnp.sqrt(2))  # Width parameters
params_init = params_init.at[0, :, 2].set(2.0)  # Set mu to (2,...,2)
for i in range(1, n):
    params_init = params_init.at[i, :, 2].set(
        2 + np.random.uniform(-0.5, 0.5, (D,))
    )  # Move to the right
# ---------- build WF + solver ---------------------------------------
oscillator = generalPotentialSolver(params_init, example_string)
SHH2 = make_SHH2_pure(splitting_type=splitting_type)  # The function to calculate S, H, H2
Smat = oscillator.calculate_S()
coeffs_init = jnp.array([1.0] * n, dtype=jnp.complex128)  # coeffs
coeffs_init = jnp.zeros_like(coeffs_init)
coeffs_init = coeffs_init.at[0].set(1.0)
norm = coeffs_init.T @ Smat @ coeffs_init
coeffs_init = coeffs_init / jnp.sqrt(norm)
initial_error = rothe_error(params_init, params_init, coeffs_init, SHH2, 0, dt)
rothesolver = RotheSolver(
    SHH2,
    dt,
    0,
    params_init,
    coeffs_init,
    rothe_grad_jax=rothe_vg_jit,
    rothe_nograd=rothe_error,
    splitting_type=splitting_type,
    name="HenonHeiles_%dD_%d" % (D, n),
    polynomial_string=example_string,
    out_dir="./wave_function_data",
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
        params=params_init,
        coeffs=coeffs_init,
    )
    print("No prior file. Starting fresh at t=0.")

rothesolver.propagate(nsteps, maxiter=1000)
# print("Using %d Gaussians in %dD Henon-Heiles potential" % (n, D))
# rothesolver.evaluate_gradient_n_times(params_init=params_init, coeffs_init=coeffs_init, n=1000)
