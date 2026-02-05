import ast
import re
import sys

import sympy as sp
from numpy import cos, exp, pi, sin

import jax.numpy as jnp


def clean_up_lines(lines: list[str]) -> list[str]:
    """
    Remove empty lines from the given list of lines.
    Also removes comments (anything after a '#').
    Also converts to lowercase.
    """
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
    for line_idx, raw in enumerate(poly_terms):  # outer index ≠ 'i'
        splitted = raw.split(":")
        poly, val = splitted
        val = val.strip()
        try:
            float(val)
            val = float(val)
        except:
            expr = sp.sympify(val)
            f = sp.lambdify(t, expr, modules="jax")
            val = f
        line_vals = []
        for d in range(dimension):  # use 'd' not 'i'
            hits = re.findall(rf"x{d}(?!\d)", poly)  # exact matches only
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
        tuple: A tuple containing the dimension (int), polynomial lines (list of str),
               and exponential lines (list of str).
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
        constants_pos = len(lines)  # Set to end if not found
    if poly_pos == -1 and exponential_pos == -1:
        raise ValueError("No 'polynomial' or 'exponential' section found in the string.")
    assert (
        constants_pos >= exponential_pos and constants_pos >= poly_pos
    ) or constants_pos == -1, "Constants section must be at the end."
    exponential_terms = []
    poly_terms = []
    constants = []
    if poly_pos < exponential_pos and exponential_pos != -1:
        # If polynomial comes first
        if poly_pos != -1:  # And if polynomial section exists
            poly_terms = lines[
                poly_pos + 1 : exponential_pos if exponential_pos != -1 else constants_pos
            ]
        exponential_terms = lines[exponential_pos + 1 : constants_pos]
    elif exponential_pos < poly_pos and poly_pos != -1:  # If exponential comes first
        poly_terms = lines[poly_pos + 1 : constants_pos]
        if exponential_pos != -1:  # If exponential section exists
            exponential_terms = lines[exponential_pos + 1 : poly_pos]
    constant_lines = lines[constants_pos + 1 :]
    poly_terms, exponential_terms = insert_constants(
        poly_terms, exponential_terms, constant_lines
    )
    return dimension, poly_terms, exponential_terms


def insert_constants(
    poly_terms: list[str], exponential_terms: list[str], constant_lines: list[str]
) -> tuple[list[str], list[str]]:
    """
    Replace constants in polynomial and exponential terms with their values.

    Args:
        poly_terms (list[str]): List of polynomial term strings.
        exponential_terms (list[str]): List of exponential term strings.
        constant_lines (list[str]): List of constant definition strings.

    Returns:
        tuple: Updated polynomial and exponential term lists.
    """
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
        raise TypeError(
            f"Exponential position must be list or tuple, got {type(pos_literal).__name__}"
        )
    pos = [float(p) for p in pos_literal]
    return lincoeff, width, pos


def make_exponential(exp_terms, dimension):
    exponential = []
    for line_idx, raw in enumerate(exp_terms):
        linCoeff, width, pos = get_exponential_coeffs(raw)
        if len(pos) != dimension:
            raise ValueError(
                f"Exponential position length does not match dimension at line {line_idx + 3}: {raw}"
            )
        exponential.append((linCoeff, width, pos))
    return exponential


def read_string(string: str):
    """
    Read a string containing polynomial coefficients and return a dictionary.

    The string should contain lines with the format:
    "x0: 1", "x1x2: 2", etc.

    Returns:
        dict: A dictionary with keys as variable names and values as coefficients.
    """
    dimension, poly_terms, exponential_terms = split_up_string(string)
    polynomial = make_polynomial(poly_terms, dimension)
    exponential = make_exponential(exponential_terms, dimension)
    return polynomial, exponential


def obtain_polynomial_squared(polynomial, t=0):
    """
    Obtain the squared polynomial coefficients from the given polynomial.

    Args:
        polynomial (list): List of lists containing variable counts and coefficients.

    Returns:
        list: List of lists with squared polynomial coefficients.
    """
    result = {}
    for i, (vars1, coeff1) in enumerate(polynomial.items()):
        for j, (vars2, coeff2) in enumerate(polynomial.items()):
            if i < j:
                continue
            duplicator = (
                1 if i == j else 2
            )  # Take into account double-counting of non-diagonal terms
            new_vars = [v1 + v2 for v1, v2 in zip(vars1, vars2)]
            c1 = coeff1 if not callable(coeff1) else float(coeff1(t))
            c2 = coeff2 if not callable(coeff2) else float(coeff2(t))

            new_coeff = c1 * c2 * duplicator
            result[tuple(new_vars)] = result.get(tuple(new_vars), 0) + new_coeff
    return result


if __name__ == "__main__":
    example_string = """
        dimension 4
        exponential
        1.0, 0.5, [0.0, 0.0, 0.0, 0.0]
        polynomial
        x0x1: 2
        x2x3: 1
        x0x2: 1
        x1x3: -1
        x0: E0*sin(pi**2/N_T*t)**2*sin(omega*t) # Example of a time-dependent string
        constants
        omega: 0.057
        E0: 0.06
        N_T: 4
        """
    polynomial, exponential = read_string(example_string)
    polynomial_squared = obtain_polynomial_squared(polynomial)
    print("Original Polynomial Coefficients:")
    print(polynomial)
    print("Exponential terms:")
    print(exponential)
    print("Squared Polynomial Coefficients:")
    print(polynomial_squared)
