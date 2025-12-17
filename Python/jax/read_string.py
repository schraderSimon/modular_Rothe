import ast
import re


def clean_up_lines(lines: list[str]) -> list[str]:
    """
    Remove empty lines from the given list of lines.
    Also removes comments (anything after a '#').
    """
    return [line.split("#")[0].strip() for line in lines if line.strip() != ""]


def find_string(lines: list[str], target: str) -> int:
    for i, line in enumerate(lines):
        if line.strip() == target:
            return i
    return -1


def make_polynomial(poly_terms, dimension):
    poly_list = []
    polynomial = {}
    for line_idx, raw in enumerate(poly_terms):  # outer index ≠ 'i'
        poly, val = raw.split(":")
        line_vals = []
        for d in range(dimension):  # use 'd' not 'i'
            hits = re.findall(rf"x{d}(?!\d)", poly)  # exact matches only
            line_vals.append(len(hits))
        if line_vals in poly_list:
            raise ValueError(f"Duplicate polynomial found at line {line_idx + 3}: {raw}")
        poly_list.append(line_vals)
        polynomial[tuple(line_vals)] = float(val.strip())
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
    if poly_pos == -1 and exponential_pos == -1:
        raise ValueError("No 'polynomial' or 'exponential' section found in the string.")
    exponential_terms = []
    poly_terms = []
    if poly_pos < exponential_pos and exponential_pos != -1:
        if poly_pos != -1:
            poly_terms = lines[
                poly_pos + 1 : exponential_pos if exponential_pos != -1 else None
            ]
        exponential_terms = lines[exponential_pos + 1 :]
    elif exponential_pos < poly_pos and poly_pos != -1:
        poly_terms = lines[poly_pos + 1 :]
        if exponential_pos != -1:
            exponential_terms = lines[exponential_pos + 1 : poly_pos]
    return dimension, poly_terms, exponential_terms


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


def obtain_polynomial_squared(polynomial):
    """
    Obtain the squared polynomial coefficients from the given polynomial.

    Args:
        polynomial (list): List of lists containing variable counts and coefficients.

    Returns:
        list: List of lists with squared polynomial coefficients.
    """
    result = {}
    iterator = list(polynomial.items())
    for i, (vars1, coeff1) in enumerate(iterator):
        for j, (vars2, coeff2) in enumerate(iterator):
            if i < j:
                continue
            duplicator = (
                1 if i == j else 2
            )  # Take into account double-counting of non-diagonal terms
            new_vars = [v1 + v2 for v1, v2 in zip(vars1, vars2)]
            new_coeff = coeff1 * coeff2 * duplicator
            result[tuple(new_vars)] = result.get(tuple(new_vars), 0) + new_coeff
    return result


if __name__ == "__main__":
    example_string = """
        dimension 4
        polynomial
        x0x1: 2
        x2x3: 1
        x0x2: 1
        x1x3: -1
        """
    polynomial, exponential = read_string(example_string)
    polynomial_squared = obtain_polynomial_squared(polynomial)
    print("Original Polynomial Coefficients:")
    print(polynomial)
    print("Exponential terms:")
    print(exponential)
    print("Squared Polynomial Coefficients:")
    print(polynomial_squared)
