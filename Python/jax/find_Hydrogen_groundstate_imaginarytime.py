import pandas as pd

import jax.numpy as jnp
from file_handling import save_rothe_state
from general_wf import generalPotentialSolver
from jax_Rothe import RotheSolver, setUpRotheErrorAndGradient_jit

df = pd.read_csv("gaussian_Coulomb/coeffs_mu100_N20.csv")
lincoeffs = df["linear"]
width = df["nonlinear"]
hydrogen_string = """
    dimension 3
    exponential
    """
for lc, w in zip(lincoeffs, width):
    hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"
initial_linear_coeffs = [
    (-0.2005325276070134 + 0j),
    (-0.42914381204851876 + 0j),
    (-0.30467705257687383 + 0j),
    (-0.13597865079192395 + 0j),
    (-0.05015841866344768 + 0j),
    (-0.01675322282470429 + 0j),
    (-0.005045114316810095 + 0j),
    (-0.0014024978113044567 + 0j),
    (-0.0003659850609822185 + 0j),
    (-9.151289205536206e-05 + 0j),
    (-2.075945675152069e-05 + 0j),
    (-3.2602513596741156e-06 + 0j),
]


n = len(
    initial_linear_coeffs
)  # use a modestly larger basis so the ground-state energy converges to -0.5
random_widths = [
    0.27718554582694666,
    0.43691440160201783,
    0.6852203154179158,
    1.0784145470226345,
    1.710821231164854,
    2.765743401093018,
    4.588621130126977,
    7.768234078253707,
    13.307147528224515,
    22.92216415164161,
    39.596193778751754,
    67.8663025891332,
]
n = len(random_widths)
params_init = jnp.zeros((n, 3, 4), dtype=jnp.float64)

for i in range(n):
    params_init = params_init.at[i, :, 0].set(random_widths[i])  # Isotropic gaussians
hydrogen_wf = generalPotentialSolver(params_init, hydrogen_string)
Smat = hydrogen_wf.calculate_S()
rothe_error, rothe_vg_jit = setUpRotheErrorAndGradient_jit(splitting_type="none")
dt = -0.001j
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
E0 = jnp.vdot(initial_linear_coeffs, H @ initial_linear_coeffs) / norm
Esq = jnp.vdot(initial_linear_coeffs, H2 @ initial_linear_coeffs) / norm
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
    name="Hydrogen_%d" % (n),
    polynomial_string=hydrogen_string,
    out_dir="./wave_function_data",
)
print(f"Initial Rothe error: {initial_error/dt**2}")
nsteps = int(abs(100 / dt))
try:
    info = rothesolver.resume_from_file(t_start=None)  # or pick a time
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

rothesolver.propagate(nsteps, maxiter=1000)
# print("Using %d Gaussians in %dD Henon-Heiles potential" % (n, D))
# rothesolver.evaluate_gradient_n_times(params_init=params_init, coeffs_init=coeffs_init, n=1000)
