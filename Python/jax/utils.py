import jax.numpy as jnp


def eigh_canonical(
    S: jnp.ndarray, H: jnp.ndarray, eps: float = 1e-12
) -> tuple[jnp.ndarray, jnp.ndarray]:
    S = _symh(S)
    H = _symh(H)
    w, U = jnp.linalg.eigh(S)

    wmax = jnp.max(jnp.real(w))
    keep = jnp.real(w) > eps * wmax

    U1, w1 = U[:, keep], w[keep]
    X = U1 * (1.0 / jnp.sqrt(w1))[None, :]

    H_orth = X.conj().T @ H @ X

    e, V = jnp.linalg.eigh(H_orth)
    C = X @ V

    return e, C


def _symh(A: jnp.ndarray) -> jnp.ndarray:  # Hermitianize (numerical hygiene)
    return 0.5 * (A + A.conj().T)


def divide_outer(a, b):
    return jnp.divide(a[:, None], b[None, :])
