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


def orthonormalize_matrices(S, H, H2, eps=1e-12):
    """
    Transform overlap and Hamiltonian matrices to orthonormal basis.

    Performs canonical orthogonalization: transforms to a basis where
    the overlap matrix is the identity.
    """
    S_sym = _symh(S)
    H_sym = _symh(H)
    H2_sym = _symh(H2)

    # Canonical orthogonalization
    w, U = jnp.linalg.eigh(S_sym)
    wmax = jnp.max(jnp.real(w))
    keep = jnp.real(w) > eps * wmax
    U1, w1 = U[:, keep], w[keep]
    X = U1 * (1.0 / jnp.sqrt(w1))[None, :]

    # Transform matrices to orthonormal basis
    H_orth = X.conj().T @ H_sym @ X
    H2_orth = X.conj().T @ H2_sym @ X

    # Diagonalize H in orthonormal basis
    eigvals, eigvecs_orth = jnp.linalg.eigh(H_orth)

    return eigvals, eigvecs_orth, H_orth, H2_orth


def calculate_variance(state_vec, H2_orth, energy):
    expectation_H2 = jnp.real(state_vec.conj() @ H2_orth @ state_vec)
    return expectation_H2 - energy**2
