import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from libraries.general_wf import *
from libraries.utils import calculate_variance, orthonormalize_matrices

df = pd.read_csv("gaussian_Coulomb/coeffs_mu=100_N=19.csv")
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

n = 20  # use a modestly larger basis so the ground-state energy converges to -0.5
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
random_widths = np.logspace(-2, 2, n, dtype=jnp.float64)
n = len(random_widths)
intial_params = jnp.zeros((n, 3, 4), dtype=jnp.float64)

for i in range(n):
    intial_params = intial_params.at[i, :, 0].set(random_widths[i])  # Isotropic gaussians
hydrogen_wf = generalPotentialSolver(intial_params, hydrogen_string)

S = hydrogen_wf.calculate_S()
S = S + 1e-8 * jnp.eye(S.shape[0], dtype=S.dtype)  # Regularization for stability
H = hydrogen_wf.calculate_H()
H2 = hydrogen_wf.calculate_H2()
result = orthonormalize_matrices(S, H, H2, return_Smin05=True)
if len(result) == 5:
    eigvals, eigvecs_orth, H_orth, H2_orth, inverse = result
else:
    eigvals, eigvecs_orth, H_orth, H2_orth = result
    inverse = None
groundstate = eigvecs_orth[:, 0]
H2_eigenvalues = jnp.linalg.eigvalsh(H2_orth)
print("H2 eigenvalues:", H2_eigenvalues[:])
print("H eigenvalues squared:", jnp.linalg.eigvalsh(H_orth)[:] ** 2)
for i in range(5):
    # Calculate variance of ith state
    state = eigvecs_orth[:, i]
    energy = eigvals[i].real
    variance = calculate_variance(state, H2_orth, energy)
    print(f"State {i}: Energy: {energy}, Variance: {variance}")
variance = calculate_variance(groundstate, H2_orth, eigvals[0].real)
print(f"Ground state energy: {eigvals[0].real}, Variance: {variance}")
# Print best linear coefficients in original basis
C = inverse @ groundstate
print(C.T @ S @ C)

print(C.T @ H @ C)
print(C.T @ H2 @ C)
print(list(np.array(C)))
assert jnp.isclose(eigvals[0].real, -0.5, atol=1e-3)
assert variance >= 0 and variance < 1e-2
