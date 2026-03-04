"""
Parser for the polynomial + exponential potential DSL.

Turns a multiline string description into:
- A polynomial dict ``{(nx0, nx1, ...): coefficient}``
- An exponential list ``[(linear_coeff, width, [positions])]``

Example input::

    dimension 2
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    exponential
    1.0, 0.5, [0.0, 0.0]
"""

import ast
import re

import sympy as sp


def clean_up_lines(lines: list[str]) -> list[str]:
    """Remove empty lines, comments, and convert to lowercase."""
    return [line.split("#")[0].strip().lower() for line in lines if line.strip() != ""]


def find_string(lines: list[str], target: str) -> int:
    for i, line in enumerate(lines):
        if line.strip() == target:
            return i
    return -1


def make_polynomial(poly_terms, dimension):
    poly_list = []
    polynomial = {}
    t = sp.symbols("t")
    for line_idx, raw in enumerate(poly_terms):
        splitted = raw.split(":")
        poly, val = splitted
        val = val.strip()
        try:
            float(val)
            val = float(val)
        except Exception:
            expr = sp.sympify(val)
            f = sp.lambdify(t, expr, modules="jax")
            val = f
        line_vals = []
        for d in range(dimension):
            hits = re.findall(rf"x{d}(?!\d)", poly)
            line_vals.append(len(hits))
        if line_vals in poly_list:
            raise ValueError(f"Duplicate polynomial found at line {line_idx + 3}: {raw}")
        poly_list.append(line_vals)
        polynomial[tuple(line_vals)] = val
    return polynomial


def split_up_string(string: str) -> tuple[int, list[str], list[str]]:
    """
    Split the input string into polynomial and exponential sections.

    Returns:
        tuple: (dimension, polynomial_lines, exponential_lines)
    """
    lines = string.splitlines()
    lines = clean_up_lines(lines)
    if not lines[0].startswith("dimension"):
        raise ValueError("First line must specify 'dimension <int>'.")
    dimension = int(lines[0].split()[1])
    poly_pos = find_string(lines, "polynomial")
    exponential_pos = find_string(lines, "exponential")
    constants_pos = find_string(lines, "constants")
    if constants_pos == -1:
        constants_pos = len(lines)
    if poly_pos == -1 and exponential_pos == -1:
        raise ValueError("No 'polynomial' or 'exponential' section found in the string.")
    assert (
        constants_pos >= exponential_pos and constants_pos >= poly_pos
    ) or constants_pos == -1, "Constants section must be at the end."
    exponential_terms = []
    poly_terms = []
    if poly_pos < exponential_pos and exponential_pos != -1:
        if poly_pos != -1:
            poly_terms = lines[poly_pos + 1 : exponential_pos if exponential_pos != -1 else constants_pos]
        exponential_terms = lines[exponential_pos + 1 : constants_pos]
    elif exponential_pos < poly_pos and poly_pos != -1:
        poly_terms = lines[poly_pos + 1 : constants_pos]
        if exponential_pos != -1:
            exponential_terms = lines[exponential_pos + 1 : poly_pos]
    constant_lines = lines[constants_pos + 1 :]
    poly_terms, exponential_terms = insert_constants(poly_terms, exponential_terms, constant_lines)
    return dimension, poly_terms, exponential_terms


def insert_constants(
    poly_terms: list[str], exponential_terms: list[str], constant_lines: list[str]
) -> tuple[list[str], list[str]]:
    """Replace constants in polynomial and exponential terms with their values."""
    exponential_terms_new = []
    poly_terms_new = []
    constant_value_dict: dict[str, str] = {}
    for i in range(len(constant_lines)):
        line = constant_lines[i]
        constant, value = line.split(":")
        constant = constant.strip()
        value = value.strip()
        constant_value_dict[constant] = value
    try:
        constant_value_dict["pi"]
    except KeyError:
        from numpy import pi

        constant_value_dict["pi"] = str(pi)
    for constant, value in constant_value_dict.items():
        for term in poly_terms:
            newterm = term.replace(constant, value)
            poly_terms_new.append(newterm)
        for term in exponential_terms:
            newterm = term.replace(constant, value)
            exponential_terms_new.append(newterm)
        exponential_terms = exponential_terms_new
        poly_terms = poly_terms_new
        exponential_terms_new = []
        poly_terms_new = []
    return poly_terms, exponential_terms


def get_exponential_coeffs(string: str) -> tuple[float, float, list[float]]:
    data = string.split(",")
    lincoeff = data[0]
    width = data[1]
    pos = ",".join(data[2:])
    lincoeff = float(lincoeff)
    width = float(width)
    try:
        pos_literal = ast.literal_eval(pos)
    except (ValueError, SyntaxError) as exc:
        raise ValueError(f"Invalid exponential position literal: {pos}") from exc
    if not isinstance(pos_literal, (list, tuple)):
        raise TypeError(f"Exponential position must be list or tuple, got {type(pos_literal).__name__}")
    pos = [float(p) for p in pos_literal]
    return lincoeff, width, pos


def make_exponential(exp_terms, dimension):
    exponential = []
    for line_idx, raw in enumerate(exp_terms):
        linCoeff, width, pos = get_exponential_coeffs(raw)
        if len(pos) != dimension:
            raise ValueError(f"Exponential position length does not match dimension at line {line_idx + 3}: {raw}")
        exponential.append((linCoeff, width, pos))
    return exponential


def read_string(string: str):
    """
    Read a potential-definition string and return (polynomial_dict, exponential_list).

    The string should define dimension, polynomial terms, and/or exponential terms.
    """
    dimension, poly_terms, exponential_terms = split_up_string(string)
    polynomial = make_polynomial(poly_terms, dimension)
    exponential = make_exponential(exponential_terms, dimension)
    return polynomial, exponential


def obtain_polynomial_squared(polynomial, t=0):
    """
    Obtain the squared polynomial coefficients from the given polynomial.

    Returns:
        dict: ``{(n0, n1, ...): coefficient}`` for the squared polynomial.
    """
    result = {}
    for i, (vars1, coeff1) in enumerate(polynomial.items()):
        for j, (vars2, coeff2) in enumerate(polynomial.items()):
            if i < j:
                continue
            duplicator = 1 if i == j else 2
            new_vars = [v1 + v2 for v1, v2 in zip(vars1, vars2)]
            c1 = coeff1 if not callable(coeff1) else float(coeff1(t))
            c2 = coeff2 if not callable(coeff2) else float(coeff2(t))
            new_coeff = c1 * c2 * duplicator
            result[tuple(new_vars)] = result.get(tuple(new_vars), 0) + new_coeff
    return result
