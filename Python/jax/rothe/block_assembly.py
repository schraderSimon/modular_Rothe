"""Shared Rothe block formulas for A, rho, and B matrix assembly.

These helpers are backend-agnostic and work with both NumPy and JAX arrays.
"""


def shh2_triplets(SHH2, t_mid, splitting_eff, ket_bra_pairs):
    """Compute multiple (S, H, H2) triplets for a sequence of (ket, bra) pairs."""
    return tuple(
        SHH2(t=t_mid, params_ket=params_ket, params_bra=params_bra, splitting_type=splitting_eff)
        for params_ket, params_bra in ket_bra_pairs
    )


def split_oldold_frozen_triplets(S_xx, H_xx, H2_xx, n_dyn):
    """Split full old-old triplets into frozen-frozen and old-to-frozen blocks."""
    frozen_frozen = (S_xx[n_dyn:, n_dyn:], H_xx[n_dyn:, n_dyn:], H2_xx[n_dyn:, n_dyn:])
    old_to_frozen = (S_xx[:, n_dyn:], H_xx[:, n_dyn:], H2_xx[:, n_dyn:])
    return frozen_frozen, old_to_frozen


def make_A_block(S, H, H2, dt_abs_sq_4, dt_imag):
    """Compute A block: S + |dt|^2/4 * H2 - Im(dt) * H."""
    return S + dt_abs_sq_4 * H2 - dt_imag * H


def make_rho_block(S, H, H2, dt_abs_sq_4, dt_real):
    """Compute rho block: S - |dt|^2/4 * H2 - i * Re(dt) * H."""
    return S - dt_abs_sq_4 * H2 - 1j * dt_real * H


def make_B_block(S, H, H2, dt_abs_sq_4, dt_imag):
    """Compute B block: S + |dt|^2/4 * H2 + Im(dt) * H."""
    return S + dt_abs_sq_4 * H2 + dt_imag * H
