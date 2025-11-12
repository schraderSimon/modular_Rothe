import numpy as np
from numpy.linalg import slogdet, det, norm
from scipy.linalg import sqrtm, logm

def make_spd(n, jitter=1, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    M = rng.standard_normal((n, n))
    S = M @ M.T  # symmetric, PSD
    S += jitter * np.eye(n)  # push eigenvalues away from zero: SPD
    return S

def make_symmetric(n, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    R = rng.standard_normal((n, n))
    return 0.5 * (R + R.T)

def det_sqrt_via_sqrtm(A):
    S = sqrtm(A)  # principal matrix square root
    return det(S)

def det_sqrt_via_trace_log(A):
    L = logm(A)  # principal matrix logarithm
    return np.exp(0.5 * np.trace(L))

def det_sqrt_via_slogdet(A):
    # slogdet returns sign (complex unit-modulus) and log|det|
    sign, logabs = slogdet(A)
    logdet = logabs + np.log(sign)  # principal complex log(det A)
    return np.exp(0.5 * logdet)

def main(n=3, seed=0, jitter=0, rtol=1e-10, atol=1e-12):
    rng = np.random.default_rng()
    ReA = make_spd(n, jitter=jitter, rng=rng)           # SPD, real symmetric
    ImA = make_symmetric(n, rng=rng)                    # real symmetric
    A = ReA.astype(np.complex128) + 1j * ImA            # complex with Re(A) SPD

    # Compute the three versions
    d1 = det_sqrt_via_sqrtm(A)
    d2 = det_sqrt_via_trace_log(A)
    d3 = det_sqrt_via_slogdet(A)

    # Print values
    print("det(sqrtm(A))                 =", d1)
    print("exp(1/2 * tr(logm(A)))       =", d2)
    print("exp(1/2 * slogdet(A).logdet) =", d3)

    # Pairwise differences
    print("\nPairwise absolute diffs:")
    print("|d1 - d2| =", abs(d1 - d2))
    print("|d1 - d3| =", abs(d1 - d3))
    print("|d2 - d3| =", abs(d2 - d3))

    # Relative checks
    def rel_close(x, y):
        denom = max(atol, abs(x), abs(y))
        return abs(x - y) <= rtol * denom + atol

    ok12 = rel_close(d1, d2)
    ok13 = rel_close(d1, d3)
    ok23 = rel_close(d2, d3)

    # Sanity: det(sqrt(A))^2 == det(A) on the same principal branch
    dA = det(A)
    sanity = rel_close(d1**2, dA)

    print("\nChecks:")
    print("d1 == d2 ?", ok12)
    print("d1 == d3 ?", ok13)
    print("d2 == d3 ?", ok23)
    print("det(sqrt(A))^2 == det(A) ?", sanity)

    if not (ok12 and ok13 and ok23 and sanity):
        raise AssertionError("Equality check failed (branch cut or conditioning issue).")

if __name__ == "__main__":
    main()
