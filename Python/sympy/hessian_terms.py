import sys
from itertools import count

from sympy import I, IndexedBase, Poly, conjugate, cse, simplify, symbols

temp = IndexedBase("temp")
x, y = symbols("x y", real=True)  # Two different dimensions are sufficient for demonstration

a_m, b_m, mu_m, p_m = symbols("a_m b_m mu_m p_m", real=True)
a_n, b_n, mu_n, p_n = symbols("a_n b_n mu_n p_n", real=True)


def a_polynomial(a, b, mu, p, x):
    return 1 / (2 * a) - 2 * a * (x - mu) ** 2


def b_polynomial(a, b, mu, p, x):
    return -I * (x - mu) ** 2


def mu_polynomial(a, b, mu, p, x):
    return 2 * (a**2 + I * b) * (x - mu) - I * p


def p_polynomial(a, b, mu, p, x):
    return I * (x - mu)


labels = ["a", "b", "mu", "p"]


def build_terms():
    a_m_k1 = a_polynomial(a_m, b_m, mu_m, p_m, x)
    a_n_k1 = a_polynomial(a_n, b_n, mu_n, p_n, x)
    a_n_k2 = a_polynomial(a_n, b_n, mu_n, p_n, y)

    b_m_k1 = b_polynomial(a_m, b_m, mu_m, p_m, x)
    b_n_k1 = b_polynomial(a_n, b_n, mu_n, p_n, x)
    b_n_k2 = b_polynomial(a_n, b_n, mu_n, p_n, y)

    mu_m_k1 = mu_polynomial(a_m, b_m, mu_m, p_m, x)
    mu_n_k1 = mu_polynomial(a_n, b_n, mu_n, p_n, x)
    mu_n_k2 = mu_polynomial(a_n, b_n, mu_n, p_n, y)

    p_m_k1 = p_polynomial(a_m, b_m, mu_m, p_m, x)
    p_n_k1 = p_polynomial(a_n, b_n, mu_n, p_n, x)
    p_n_k2 = p_polynomial(a_n, b_n, mu_n, p_n, y)

    terms_m = [a_m_k1, b_m_k1, mu_m_k1, p_m_k1]
    terms_n_same = [a_n_k1, b_n_k1, mu_n_k1, p_n_k1]
    terms_n_diff = [a_n_k2, b_n_k2, mu_n_k2, p_n_k2]
    return terms_m, terms_n_same, terms_n_diff


labels_m = ["a_m", "b_m", "mu_m", "p_m"]
labels_n = ["a_n", "b_n", "mu_n", "p_n"]


def get_polynomial_coeffs_diff():
    terms_m, terms_n_same, terms_n_diff = build_terms()
    m_terms = []
    for m, term_m in enumerate(terms_m):
        m_term = Poly(term_m, x)
        m_terms.append(m_term)
    n_terms = []
    m_out = []
    for n, term_n in enumerate(terms_n_diff):
        n_term = Poly(conjugate(term_n), y)
        n_terms.append(n_term)
    print("M-terms:")

    for m, term_m in enumerate(terms_m):
        typelist = []
        for m_power, coeff in sorted(m_terms[m].as_dict().items()):

            # print(f"M-term {labels_m[m]}, x^{m_power[0]}: coeff = {coeff}")
            typelist.append(coeff)
        m_out.append(typelist)
    print(m_out)
    n_out = []
    print("N-terms:")
    for n, term_n in enumerate(terms_n_diff):
        typelist = []
        for n_power, coeff in sorted(n_terms[n].as_dict().items()):
            # print(f"N-term {labels_n[n]}, y^{n_power[0]}: coeff = {coeff}")
            typelist.append(coeff)
        n_out.append(typelist)
    print(n_out)


def global_cse_on_all_coeffs():
    terms_m, terms_n_same, terms_n_diff = build_terms()

    all_exprs = []  # all coefficient expressions to feed into cse
    meta = []  # bookkeeping so we can map back

    # Same-dimension terms: polynomials in x
    for m, term_m in enumerate(terms_m):
        for n, term_n in enumerate(terms_n_same):
            expr = conjugate(term_m) * term_n
            P = Poly(expr, x)
            for (power,), c in sorted(P.as_dict().items()):
                all_exprs.append(simplify(c))
                meta.append(("same", m, n, power))

    # Cross-dimension terms: polynomials in (x, y)
    """
    for m, term_m in enumerate(terms_m):
        for n, term_n in enumerate(terms_n_diff):
            expr_xy = conjugate(term_m) * term_n
            Pxy = Poly(expr_xy, x, y)
            for (ix, iy), c in sorted(Pxy.as_dict().items()):
                all_exprs.append(simplify(c))
                meta.append(("diff", m, n, ix, iy))
    """
    # One global CSE over everything
    temps, reduced = cse(
        all_exprs,
        symbols=(temp[i] for i in count()),  # temp[0], temp[1], ...
        optimizations="basic",
    )

    # Rebuild structured coefficient dictionaries
    coeffs_same = {}  # (m, n) -> {power: expr}
    coeffs_diff = {}  # (m, n) -> {(ix, iy): expr}

    for info, c_red in zip(meta, reduced):
        kind = info[0]
        if kind == "same":
            _, m, n, power = info
            key = (m, n)
            if key not in coeffs_same:
                coeffs_same[key] = {}
            coeffs_same[key][power] = c_red
        else:
            _, m, n, ix, iy = info
            key = (m, n)
            if key not in coeffs_diff:
                coeffs_diff[key] = {}
            coeffs_diff[key][(ix, iy)] = c_red

    return temps, coeffs_same, coeffs_diff


def print_results():
    temps, coeffs_same, coeffs_diff = global_cse_on_all_coeffs()

    print("Global intermediates (reused across all 32 products):")
    for s, e in temps:
        print(f"{s} = {e}")
    print()
    print("----------------------------------------")
    print("Same-dimension coefficients (polynomials in x):")
    powers = [0, 1, 2, 3, 4]
    for power in powers:
        print()
        print()
        print("Coefficients for x^", power)
        for (m, n), d in coeffs_same.items():

            try:
                d[power]
                # print(labels[m], labels[n])
                print(f"M[:,:,:,{m},{n}]={d[power]}")
            except KeyError:
                continue
    print()
    print("----------------------------------------")


if __name__ == "__main__":
    # print_results()
    get_polynomial_coeffs_diff()
