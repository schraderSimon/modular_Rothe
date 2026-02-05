import pandas as pd

import jax.numpy as jnp
from libraries.file_handling import save_rothe_state
from libraries.general_wf import generalPotentialSolver
from libraries.jax_Rothe import RotheSolver, setUpRotheErrorAndGradient_jit
from read_best_coefficients import read_best_coefficients

num_gauss_potential = 21
num_gauss_wavefunction = 50
df = pd.read_csv("gaussian_Coulomb/coeffs_mu=100_N=%d.csv" % (num_gauss_potential))
lincoeffs = df["linear"]
width = df["nonlinear"]
hydrogen_string = """
    dimension 3
    polynomial
    x2: E_0*sin(omega*t)*sin(omega*t/N_T/2)**2
    exponential
    """
#

for lc, w in zip(lincoeffs, width):
    hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"
hydrogen_string += f"constants\nN_T: 3\nomega: 0.057\nE_0: 0.06"

n = num_gauss_wavefunction
initial_linear_coeffs, params_init = read_best_coefficients(
    num_gauss_wavefunction=num_gauss_wavefunction, num_gauss_potential=num_gauss_potential
)
print(initial_linear_coeffs.shape)
print(params_init.shape)
hydrogen_wf = generalPotentialSolver(params_init, hydrogen_string)
Smat = hydrogen_wf.calculate_S()
dt = 0.2
rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type="none")

D = 3


def make_SHH2_pure(splitting_type="none"):
    """
    Wrap `generalPotentialSolver.calculate_SHH2` for Rothe.

    `splitting_type` should be "none" to include both kinetic and potential.
    Passing any other value returns only the potential part, which would miss
    the kinetic energy and give the wrong total energy.
    """

    system = generalPotentialSolver(params_init, hydrogen_string)

    def SHH2(t, params, splitting_type=splitting_type):
        system.update_parameters(params)
        return system.calculate_SHH2(t=t, splitting_type=splitting_type)

    return SHH2


# Use the full Hamiltonian (kinetic + potential).
SHH2 = make_SHH2_pure(splitting_type="none")

initial_linear_coeffs = jnp.array(initial_linear_coeffs)
S, H, H2 = SHH2(0, params_init)
norm = jnp.vdot(initial_linear_coeffs, S @ initial_linear_coeffs)
initial_linear_coeffs = initial_linear_coeffs / jnp.sqrt(norm)
E0 = jnp.vdot(initial_linear_coeffs, H @ initial_linear_coeffs)
Esq = jnp.vdot(initial_linear_coeffs, H2 @ initial_linear_coeffs)
print("Initial energy:", E0)
print("Initial norm:", norm)
print("Initial variance:", Esq - E0**2)
initial_error = rothe_error(params_init, params_init, initial_linear_coeffs, SHH2, 0, dt)
rothesolver = RotheSolver(
    SHH2,
    dt,
    0,
    params_init,
    initial_linear_coeffs,
    rothe_grad_jax=rothe_vg_jit,
    rothe_nograd=rothe_error,
    splitting_type="none",
    name="Hydrogen_%d_%d" % (n, num_gauss_potential),
    polynomial_string=hydrogen_string,
    out_dir="./wave_function_data",
)
nsteps = int(abs(1000 / dt))
try:
    info = rothesolver.resume_from_file(t_start=None)  # resume from latest saved time
    print(f"Resumed at t={info['t']} (step {info['idx']}), trimmed={info['trimmed']}")
    nsteps = nsteps - info["idx"]  # remaining steps
except FileNotFoundError:
    # optional: log initial state at t=0, so resumes have an anchor
    save_rothe_state(
        rothesolver.name,
        "none",
        hydrogen_string,
        dt,
        0.0,
        float(initial_error),
        params=params_init,
        coeffs=initial_linear_coeffs,
    )
    print("No prior file. Starting fresh at t=0.")

rothesolver.propagate(nsteps, maxiter=300)
# print("Using %d Gaussians in %dD Henon-Heiles potential" % (n, D))
# rothesolver.evaluate_gradient_n_times(params_init=params_init, coeffs_init=coeffs_init, n=1000)
