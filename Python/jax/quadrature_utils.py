import numpy as np


def eval_gaussian_basis(params, idx, X, Y):
    """Evaluate normalized basis Gaussian g_idx(x) on a 2D grid (supports b,p)."""
    a = np.asarray(params[idx, :, 0])
    b = np.asarray(params[idx, :, 1])
    mu = np.asarray(params[idx, :, 2])
    p = np.asarray(params[idx, :, 3])
    pref = np.prod((2 * a**2 / np.pi) ** 0.25)
    expo = -(a[0] ** 2 + 1j * b[0]) * (X - mu[0]) ** 2
    expo -= (a[1] ** 2 + 1j * b[1]) * (Y - mu[1]) ** 2
    expo += 1j * (p[0] * (X - mu[0]) + p[1] * (Y - mu[1]))
    return pref * np.exp(expo)


def eval_gaussian_potential(exp_params, lin_params, X, Y):
    """Evaluate summed Gaussian potential on grid (supports b,p)."""
    total = np.zeros_like(X, dtype=np.complex128)
    for k in range(exp_params.shape[0]):
        a = np.asarray(exp_params[k, :, 0])
        b = np.asarray(exp_params[k, :, 1])
        mu = np.asarray(exp_params[k, :, 2])
        p = np.asarray(exp_params[k, :, 3])
        c = np.asarray(lin_params[k])
        expo = -(a[0] ** 2 + 1j * b[0]) * (X - mu[0]) ** 2
        expo -= (a[1] ** 2 + 1j * b[1]) * (Y - mu[1]) ** 2
        expo += 1j * (p[0] * (X - mu[0]) + p[1] * (Y - mu[1]))
        total += c * np.exp(expo)
    return total


def eval_polynomial(poly_dict, X, Y):
    """Evaluate polynomial defined by monomial-power keys on grid."""
    res = np.zeros_like(X, dtype=np.complex128)
    for key, coeff in poly_dict.items():
        term = coeff
        if key[0] > 0:
            term = term * (X ** key[0])
        if key[1] > 0:
            term = term * (Y ** key[1])
        res += term
    return res


def quadrature_cross_terms(params, poly_dict, exp_params, lin_params, xs, ys):
    """Compute 2*<g_i|V_gauss V_poly|g_j> via 2D quadrature."""
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    X, Y = np.meshgrid(xs, ys, indexing="ij")

    V_poly = eval_polynomial(poly_dict, X, Y)
    V_gauss = eval_gaussian_potential(exp_params, lin_params, X, Y)

    n = params.shape[0]
    numeric = np.zeros((n, n), dtype=np.complex128)
    for i in range(n):
        gi = eval_gaussian_basis(params, i, X, Y)
        for j in range(n):
            gj = eval_gaussian_basis(params, j, X, Y)
            integrand = np.conj(gi) * V_gauss * V_poly * gj
            numeric[i, j] = 2.0 * np.sum(integrand) * dx * dy
    return numeric


def eval_basis_grad_lap(params, idx, X, Y):
    a = np.asarray(params[idx, :, 0])
    b = np.asarray(params[idx, :, 1])
    mu = np.asarray(params[idx, :, 2])
    p = np.asarray(params[idx, :, 3])
    g = eval_gaussian_basis(params, idx, X, Y)

    yx, yy = X - mu[0], Y - mu[1]
    alpha0 = a[0] ** 2 + 1j * b[0]
    alpha1 = a[1] ** 2 + 1j * b[1]

    fx = -2 * alpha0 * yx + 1j * p[0]
    fy = -2 * alpha1 * yy + 1j * p[1]
    gx = fx * g
    gy = fy * g

    lapx = (-2 * alpha0 + fx**2) * g
    lapy = (-2 * alpha1 + fy**2) * g
    lap = lapx + lapy
    return g, (gx, gy), lap


def eval_potential_grad_lap(exp_params, lin_params, X, Y):
    V = np.zeros_like(X, dtype=np.complex128)
    Gx = np.zeros_like(X, dtype=np.complex128)
    Gy = np.zeros_like(X, dtype=np.complex128)
    Lap = np.zeros_like(X, dtype=np.complex128)

    for k in range(exp_params.shape[0]):
        a = np.asarray(exp_params[k, :, 0])
        mu = np.asarray(exp_params[k, :, 2])
        c = np.asarray(lin_params[k])

        yx, yy = X - mu[0], Y - mu[1]
        alpha0 = a[0] ** 2
        alpha1 = a[1] ** 2
        Vk = c * np.exp(-alpha0 * yx**2 - alpha1 * yy**2)

        gx = -2 * alpha0 * yx * Vk
        gy = -2 * alpha1 * yy * Vk
        lapk = (
            (4 * alpha0**2 * yx**2 - 2 * alpha0) + (4 * alpha1**2 * yy**2 - 2 * alpha1)
        ) * Vk

        V += Vk
        Gx += gx
        Gy += gy
        Lap += lapk

    return V, (Gx, Gy), Lap


def quadrature_gaussian_kinetic_cross(params, exp_params, lin_params, xs, ys):
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    X, Y = np.meshgrid(xs, ys, indexing="ij")

    V, (Gx, Gy), LapV = eval_potential_grad_lap(exp_params, lin_params, X, Y)

    n = params.shape[0]
    numeric = np.zeros((n, n), dtype=np.complex128)
    for i in range(n):
        gi, (gix, giy), lapi = eval_basis_grad_lap(params, i, X, Y)
        for j in range(n):
            gj, (gjx, gjy), lapj = eval_basis_grad_lap(params, j, X, Y)
            grad_term = Gx * gjx + Gy * gjy
            integrand = np.conj(gi) * (-0.5 * LapV * gj - grad_term - V * lapj)
            numeric[i, j] = np.sum(integrand) * dx * dy
    # Return TV + VT = TV + TV^†
    return numeric + numeric.conj().T
