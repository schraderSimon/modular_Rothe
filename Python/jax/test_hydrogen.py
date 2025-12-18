import pandas as pd

from general_wf import *
from utils import calculate_variance, orthonormalize_matrices

df = pd.read_csv("gaussian_Coulomb/coeffs_mu100_N15.csv")
lincoeffs = df["linear"]
width = df["nonlinear"]
x = np.linspace(0, 10, 1000)


def lincombGauss(linear, width, x):
    # This is the correct way to interpret the linear combination of Gaussians
    def gauss(lin, w):
        return lin * jnp.exp(-(w**2) * x**2)

    total = jnp.zeros_like(x)
    for l, w in zip(linear, width):
        total += gauss(l, w)
    return total


hydrogen_string = """
    dimension 3
    exponential
    """
for lc, w in zip(lincoeffs, width):
    hydrogen_string += f"{-lc}, {w}, [0.0, 0.0, 0.0]\n"

n = 30  # use a modestly larger basis so the ground-state energy converges to -0.5
random_widths = np.logspace(-3, 3, n)
intial_params = jnp.zeros((n, 3, 4), dtype=jnp.float64)

for i in range(n):
    intial_params = intial_params.at[i, :, 0].set(random_widths[i])  # Isotropic gaussians
hydrogen_wf = generalPotentialSolver(intial_params, hydrogen_string)

S = hydrogen_wf.calculate_S()
H = hydrogen_wf.calculate_H()
H2 = hydrogen_wf.calculate_H2()
eigvals, eigvecs_orth, H_orth, H2_orth = orthonormalize_matrices(S, H, H2)
groundstate = eigvecs_orth[:, 0]
variance = calculate_variance(groundstate, H2_orth, eigvals[0].real)
print(f"Ground state energy: {eigvals[0].real}, Variance: {variance}")
assert jnp.isclose(eigvals[0].real, -0.5, atol=1e-3)
assert variance >= 0 and variance < 1e-2
