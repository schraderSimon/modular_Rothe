import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg
import scipy.optimize
from scipy.special import erf

PROJECT_ROOT = Path(__file__).parent.parent


def get_Vgauss(lc, a):
    def gauss_func(r):
        return np.sum(-lc[:, None] * np.exp(-(a[:, None] ** 2) * r[None, :] ** 2), axis=0)

    return gauss_func


def plot(gauss_func):
    r = np.logspace(-5, 3, 1000)
    V_gauss = gauss_func(r)
    V_gauss_squared = V_gauss**2
    V_exact = erf(100 * r) ** 2 / r**2
    plt.plot(r, V_gauss_squared, label="Gaussian fit (squared)")
    plt.plot(r, V_exact, "--", label=r"$(erf(100r)/r)**2$")
    plt.plot(r, abs(V_gauss_squared - V_exact) / V_exact, label="Relative error vs exact")
    plt.xlabel("r")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel("V^2(r)")
    plt.legend()
    plt.tight_layout()
    plt.show()


def predict_yhat(r_grid, beta, w):
    r2 = r_grid**2
    A = np.exp(-r2[:, None] * beta[None, :])
    return A @ w


def solve_linear_weights(beta, r_grid, y, ridge=1e-14, weight_mode="rel", rel_floor=1e-12):
    """
    Variable projection inner solve:
        w*(beta) = argmin_w || W (A(beta) w - y) ||_2^2 + ridge ||w||_2^2

    IMPORTANT: for relative-error targeting we want residuals scaled by 1/|y|,
    i.e. W_kk = 1/max(|y_k|, eps). That corresponds to wgt = 1/denom^2 here.
    """
    beta = np.asarray(beta, dtype=float)
    beta = np.maximum(beta, 1e-300)

    r2 = r_grid**2
    A = np.exp(-r2[:, None] * beta[None, :])

    if weight_mode == "abs":
        y_w = y
        A_w = A
        sw = None
    else:
        ymax = np.max(np.abs(y))
        if weight_mode == "rel":
            denom = np.maximum(np.abs(y), 1e-300)
        elif weight_mode == "hybrid":
            denom = np.maximum(np.abs(y), rel_floor * ymax)
            denom = np.maximum(denom, 1e-300)
        else:
            raise ValueError(f"Unknown weight_mode={weight_mode}")

        # Make objective approximate sum( |err|^2 / denom^2 )  -> squared relative error
        wgt = 1.0 / (denom**2)
        sw = np.sqrt(wgt)  # = 1/denom
        y_w = y * sw
        A_w = A * sw[:, None]

    # Column scaling + ridge for stability
    col_norms = np.linalg.norm(A_w, axis=0)
    col_norms = np.maximum(col_norms, 1e-300)
    A_sc = A_w / col_norms[None, :]

    ATA = A_sc.T @ A_sc
    ATy = A_sc.T @ y_w
    w_sc = np.linalg.solve(ATA + ridge * np.eye(ATA.shape[0]), ATy)
    w = w_sc / col_norms

    y_hat = A @ w
    r = y_hat - y

    if weight_mode == "abs":
        obj = float(np.dot(r, r))
    else:
        rr = r * sw  # (err / denom)
        obj = float(np.dot(rr, rr))
    return w, y_hat, obj


def compress_square_from_gauss_sum(a, c, r_grid, M=200, tol=1e-16):
    """
    Compress g(r)^2 where g(r)=sum_i c_i exp(-(a_i^2) r^2) into M Gaussians:
        g(r)^2 ≈ sum_m w_m exp(-beta_m r^2)

    by selecting a subset of exact pairwise-sum rates beta_ij=a_i^2+a_j^2 via pivoted QR.
    """
    a2 = a**2
    beta_full = (a2[:, None] + a2[None, :]).reshape(-1)  # (N^2,)
    w_full = (c[:, None] * c[None, :]).reshape(-1)  # (N^2,)

    keep = np.abs(w_full) > tol * np.max(np.abs(w_full))
    beta_full = beta_full[keep]
    w_full = w_full[keep]

    r2 = r_grid**2
    A = np.exp(-r2[:, None] * beta_full[None, :])
    y = A @ w_full

    Q, R, piv = scipy.linalg.qr(A, pivoting=True, mode="economic")
    sel = piv[:M]
    beta_sel = beta_full[sel]

    A_sel = A[:, sel]
    w_sel, *_ = np.linalg.lstsq(A_sel, y, rcond=None)

    return np.sort(beta_sel), w_sel


def error_metrics(y_hat, y):
    denom = np.maximum(np.abs(y), 1e-300)
    rel = np.abs(y_hat - y) / denom
    max_rel = float(np.max(rel))
    rms_rel = float(np.sqrt(np.mean(rel**2)))

    scale = float(np.max(np.abs(y)))
    mask = np.abs(y) > 1e-12 * scale
    if np.any(mask):
        rel_m = np.abs(y_hat[mask] - y[mask]) / np.abs(y[mask])
        max_rel_m = float(np.max(rel_m))
        rms_rel_m = float(np.sqrt(np.mean(rel_m**2)))
    else:
        max_rel_m = float("nan")
        rms_rel_m = float("nan")

    return max_rel, rms_rel, max_rel_m, rms_rel_m


def beta_bounds_from_grid(r_grid, beta0):
    rmin2 = float(np.min(r_grid) ** 2)
    rmax2 = float(np.max(r_grid) ** 2)

    beta_lo_grid = 1e-3 / rmax2
    beta_hi_grid = 60.0 / rmin2  # a bit looser than before

    beta0 = np.asarray(beta0, dtype=float)
    beta0 = np.maximum(beta0, 1e-300)

    beta_lo = max(beta_lo_grid, float(np.min(beta0)) / 1000.0, 1e-300)
    beta_hi = min(beta_hi_grid, float(np.max(beta0)) * 1000.0)
    beta_hi = max(beta_hi, beta_lo * 10.0)
    return beta_lo, beta_hi


def optimize_betas_coord_desc(
    r_grid_opt, y_opt, beta0, ridge=1e-12, weight_mode="rel", rel_floor=1e-12, sweeps=2, min_log_gap=1e-4
):
    """
    Derivative-free coordinate descent on log(beta) with 1D bounded minimization.
    Keeps betas sorted and enforces a minimum spacing in log-space to prevent collapse.
    """
    beta0 = np.sort(np.maximum(np.asarray(beta0, dtype=float), 1e-300))
    beta_lo, beta_hi = beta_bounds_from_grid(r_grid_opt, beta0)

    logb = np.log(beta0)
    global_lo = np.log(beta_lo)
    global_hi = np.log(beta_hi)

    # baseline
    w, y_hat, best_obj = solve_linear_weights(
        beta=np.exp(logb), r_grid=r_grid_opt, y=y_opt, ridge=ridge, weight_mode=weight_mode, rel_floor=rel_floor
    )

    for _ in range(sweeps):
        for m in range(len(logb)):
            lo = global_lo
            hi = global_hi
            if m > 0:
                lo = max(lo, logb[m - 1] + min_log_gap)
            if m < len(logb) - 1:
                hi = min(hi, logb[m + 1] - min_log_gap)
            if not (lo < hi):
                continue

            def obj_1d(x):
                tmp = logb.copy()
                tmp[m] = x
                beta = np.exp(tmp)
                _, _, obj = solve_linear_weights(
                    beta=beta, r_grid=r_grid_opt, y=y_opt, ridge=ridge, weight_mode=weight_mode, rel_floor=rel_floor
                )
                return obj

            res = scipy.optimize.minimize_scalar(obj_1d, bounds=(lo, hi), method="bounded")
            if np.isfinite(res.fun) and res.fun < best_obj:
                logb[m] = float(res.x)
                best_obj = float(res.fun)

    beta_opt = np.exp(logb)
    w_opt, y_hat_opt, obj_opt = solve_linear_weights(
        beta=beta_opt, r_grid=r_grid_opt, y=y_opt, ridge=ridge, weight_mode=weight_mode, rel_floor=rel_floor
    )
    return beta_opt, w_opt, y_hat_opt, obj_opt


def main(plot_func=False):
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=31, help="Number of Gaussians in V_gauss")
    parser.add_argument("--M", type=int, default=80, help="Number of Gaussians for the squared fit")
    parser.add_argument("--grid_lo", type=float, default=-6.0, help="log10(r_min) for full grid")
    parser.add_argument("--grid_hi", type=float, default=3.0, help="log10(r_max) for full grid")
    parser.add_argument("--K", type=int, default=2000, help="Number of points on full grid")

    parser.add_argument("--K_opt", type=int, default=400, help="Number of points on optimization grid")
    parser.add_argument("--tol_drop", type=float, default=1e-16, help="Drop threshold for tiny pair weights")

    parser.add_argument("--ridge", type=float, default=1e-12, help="Ridge for linear solve in varpro")
    parser.add_argument("--opt", action="store_true", help="Refine beta via coordinate descent (non-gradient)")
    parser.add_argument("--sweeps", type=int, default=2, help="Coordinate-descent sweeps over all betas")
    parser.add_argument("--min_log_gap", type=float, default=1e-4, help="Min spacing between log betas")

    parser.add_argument(
        "--weight_mode",
        type=str,
        default="rel",
        choices=["abs", "rel", "hybrid"],
        help="Weighting for varpro objective",
    )
    parser.add_argument("--rel_floor", type=float, default=1e-12, help="Floor for hybrid weighting")
    parser.add_argument("--plot", action="store_true", help="Show diagnostic plots")
    args = parser.parse_args()

    csv_path = PROJECT_ROOT / "gaussian_Coulomb" / f"coeffs_mu=100_N={args.n}.csv"
    df = pd.read_csv(csv_path)
    lc, a = df["linear"].values, df["nonlinear"].values
    a = np.abs(a)

    func = get_Vgauss(lc, a)
    if args.plot:
        plot(func)

    grid = np.logspace(args.grid_lo, args.grid_hi, args.K)
    V_gauss = func(grid)
    V_gauss_squared = V_gauss**2

    # smaller grid used only for nonlinear optimization (much faster and usually enough)
    grid_opt = np.logspace(args.grid_lo, args.grid_hi, args.K_opt)
    V_gauss_opt = func(grid_opt)
    V_gauss_sq_opt = V_gauss_opt**2

    print(float(a.min()))
    print(float(a.max()))

    beta_init, w_init = compress_square_from_gauss_sum(a, lc, grid_opt, M=args.M, tol=args.tol_drop)

    # Evaluate initial fit on FULL grid (weights re-solved there, because why not)
    w_full_init, yhat_full_init, obj_full_init = solve_linear_weights(
        beta=beta_init,
        r_grid=grid,
        y=V_gauss_squared,
        ridge=args.ridge,
        weight_mode=args.weight_mode,
        rel_floor=args.rel_floor,
    )
    max_rel, rms_rel, max_rel_m, rms_rel_m = error_metrics(yhat_full_init, V_gauss_squared)
    print("---- Init (subset betas + varpro weights) ----")
    print("n_beta =", len(beta_init))
    print("objective(full) =", obj_full_init)
    print("max relative error =", max_rel)
    print("rms relative error =", rms_rel)
    print("masked max rel err =", max_rel_m)
    print("masked rms rel err =", rms_rel_m)

    beta_best, w_best, yhat_best = beta_init, w_full_init, yhat_full_init

    if args.opt:
        beta_opt, w_opt, yhat_opt, obj_opt = optimize_betas_coord_desc(
            r_grid_opt=grid_opt,
            y_opt=V_gauss_sq_opt,
            beta0=beta_init,
            ridge=args.ridge,
            weight_mode=args.weight_mode,
            rel_floor=args.rel_floor,
            sweeps=args.sweeps,
            min_log_gap=args.min_log_gap,
        )

        # After optimizing betas on opt-grid, re-solve weights on FULL grid
        w_full, yhat_full, obj_full = solve_linear_weights(
            beta=beta_opt,
            r_grid=grid,
            y=V_gauss_squared,
            ridge=args.ridge,
            weight_mode=args.weight_mode,
            rel_floor=args.rel_floor,
        )
        beta_best, w_best, yhat_best = beta_opt, w_full, yhat_full

        max_rel, rms_rel, max_rel_m, rms_rel_m = error_metrics(yhat_best, V_gauss_squared)
        print("---- After beta optimization (coord desc) ----")
        print("objective(opt-grid) =", obj_opt)
        print("objective(full)     =", obj_full)
        print("n_beta =", len(beta_best))
        print("max relative error =", max_rel)
        print("rms relative error =", rms_rel)
        print("masked max rel err =", max_rel_m)
        print("masked rms rel err =", rms_rel_m)

    # Write results to CSV in the same directory as the source coefficients.
    # "nonlinear" column stores sqrt(beta) so that the convention matches the
    # original files (where the stored value is squared to form the exponent).
    out_path = PROJECT_ROOT / "gaussian_Coulomb" / f"square_coeffs_mu=100_Norig={args.n}_Nfit={len(beta_best)}.csv"
    pd.DataFrame({"linear": w_best, "nonlinear": np.sqrt(beta_best)}).to_csv(out_path, index=False)
    print(f"Wrote {out_path}")

    if args.plot:
        plt.figure(figsize=(8, 5))
        plt.plot(grid, V_gauss_squared, label="target: V_gauss^2")
        plt.plot(grid, yhat_best, "--", label="fit: sum w exp(-beta r^2)")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("r")
        plt.ylabel("V^2(r)")
        plt.legend()
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(8, 5))
        denom = np.maximum(np.abs(V_gauss_squared), 1e-300)
        rel = np.abs(yhat_best - V_gauss_squared) / denom
        plt.plot(grid, rel, label="pointwise relative error")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("r")
        plt.ylabel("rel err")
        plt.legend()
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
