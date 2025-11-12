import jax.numpy as jnp
from HOscillator_1d import HOscillator, groundstate_energy
from utils import eigh_canonical


def test_implementation_1D():
    params = jnp.array([[1 / jnp.sqrt(2), 0.0, 0.0, 0.0]])  # Ground state of HO
    hoscillator = HOscillator(params)
    S, H, H2 = hoscillator.calculate_SHH2()
    energy = H[0, 0]
    assert jnp.isclose(energy, 0.5)
    variance = H2[0, 0] - energy**2
    assert jnp.isclose(variance, 0.0)
    assert jnp.isclose(groundstate_energy(params), 0.5)
    num = 20
    params = jnp.zeros((num, 4))
    params = params.at[:, 0].set(1 / jnp.sqrt(10))  # Width parameters
    params = params.at[: num // 2, 3].set(
        jnp.linspace(-3, 3, num // 2)
    )  # First dimension gets different mu parameter
    params = params.at[: num // 2, 4].set(0.1)  # Shift in momentum
    params = params.at[: num // 2, 1].set(0.1)
    params = params.at[num // 2 :, 3].set(jnp.linspace(3, -3, num // 2))
    params = params.at[num // 2 :, 4].set(-0.1)
    params = params.at[num // 2 :, 1].set(-0.1)
    hoscillator = HOscillator(params)
    S, H, H2 = hoscillator.calculate_SHH2()
    eigvals, eigvecs = eigh_canonical(S, H)
    print(eigvals)
    for i in range(3):
        assert jnp.isclose(eigvals[i], i + 0.5, rtol=1e-2)


test_implementation_1D()
