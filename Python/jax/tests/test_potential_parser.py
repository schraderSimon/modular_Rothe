"""Tests for rothe.potential_parser."""

import pytest

from rothe.potential_parser import (
    clean_up_lines,
    find_string,
    get_exponential_coeffs,
    insert_constants,
    make_exponential,
    make_polynomial,
    obtain_polynomial_squared,
    read_string,
    split_up_string,
)

# ---------------------------------------------------------------------------
# clean_up_lines
# ---------------------------------------------------------------------------


def test_clean_up_lines_removes_empty_lines():
    lines = ["", "dimension 2", "  ", "x0x0: 0.5"]
    assert clean_up_lines(lines) == ["dimension 2", "x0x0: 0.5"]


def test_clean_up_lines_strips_inline_comments():
    lines = ["polynomial  # heading", "x0x0: 0.5  # kinetic term"]
    assert clean_up_lines(lines) == ["polynomial", "x0x0: 0.5"]


def test_clean_up_lines_lowercases():
    lines = ["Dimension 2", "Polynomial"]
    assert clean_up_lines(lines) == ["dimension 2", "polynomial"]


# ---------------------------------------------------------------------------
# find_string
# ---------------------------------------------------------------------------


def test_find_string_returns_correct_index():
    lines = ["dimension 1", "polynomial", "x0x0: 0.5"]
    assert find_string(lines, "polynomial") == 1


def test_find_string_returns_minus_one_when_absent():
    lines = ["dimension 1", "x0x0: 0.5"]
    assert find_string(lines, "exponential") == -1


# ---------------------------------------------------------------------------
# split_up_string
# ---------------------------------------------------------------------------


def test_split_up_string_parses_dimension():
    s = "dimension 3\npolynomial\nx0x0: 1.0\n"
    dim, _, _ = split_up_string(s)
    assert dim == 3


def test_split_up_string_requires_dimension_as_first_line():
    with pytest.raises(ValueError, match="First line must specify"):
        split_up_string("polynomial\nx0x0: 0.5\n")


def test_split_up_string_polynomial_only_yields_no_exp_terms():
    s = "dimension 2\npolynomial\nx0x0: 0.5\nx1x1: 0.5\n"
    _, poly_terms, exp_terms = split_up_string(s)
    assert exp_terms == []
    assert len(poly_terms) == 2


def test_split_up_string_exponential_only_yields_no_poly_terms():
    s = "dimension 1\nexponential\n1.0, 0.5, [0.0]\n"
    _, poly_terms, exp_terms = split_up_string(s)
    assert poly_terms == []
    assert len(exp_terms) == 1


# ---------------------------------------------------------------------------
# make_polynomial
# ---------------------------------------------------------------------------


def test_make_polynomial_1d_quadratic():
    poly = make_polynomial(["x0x0: 0.5"], dimension=1)
    assert poly == {(2,): 0.5}


def test_make_polynomial_2d_diagonal():
    poly = make_polynomial(["x0x0: 0.5", "x1x1: 0.5"], dimension=2)
    assert poly == {(2, 0): 0.5, (0, 2): 0.5}


def test_make_polynomial_cross_term_2d():
    # "x0x1" contributes one power of x0 and one of x1 → key (1, 1)
    poly = make_polynomial(["x0x1: 1.0"], dimension=2)
    assert poly == {(1, 1): 1.0}


def test_make_polynomial_cubic_term():
    # "x0x0x0" → three factors of x0 → key (3,)
    poly = make_polynomial(["x0x0x0: 2.0"], dimension=1)
    assert poly == {(3,): 2.0}


def test_make_polynomial_duplicate_term_raises():
    with pytest.raises(ValueError, match="Duplicate polynomial"):
        make_polynomial(["x0x0: 0.5", "x0x0: 1.0"], dimension=1)


# ---------------------------------------------------------------------------
# get_exponential_coeffs
# ---------------------------------------------------------------------------


def test_get_exponential_coeffs_parses_correctly():
    lincoeff, width, pos = get_exponential_coeffs("1.5, 2.0, [0.5, -0.5]")
    assert lincoeff == pytest.approx(1.5)
    assert width == pytest.approx(2.0)
    assert pos == pytest.approx([0.5, -0.5])


def test_get_exponential_coeffs_set_position_raises_for_non_list():
    # A dict literal is not a valid list/tuple position
    with pytest.raises(TypeError, match="list or tuple"):
        get_exponential_coeffs("1.0, 0.5, {0.0: 1}")


# ---------------------------------------------------------------------------
# make_exponential
# ---------------------------------------------------------------------------


def test_make_exponential_correct_output():
    result = make_exponential(["1.0, 0.5, [0.0, 0.0]"], dimension=2)
    assert len(result) == 1
    lincoeff, width, pos = result[0]
    assert lincoeff == pytest.approx(1.0)
    assert width == pytest.approx(0.5)
    assert pos == pytest.approx([0.0, 0.0])


def test_make_exponential_dimension_mismatch_raises():
    # 2D position vector but dimension=1
    with pytest.raises(ValueError, match="does not match dimension"):
        make_exponential(["1.0, 0.5, [0.0, 0.0]"], dimension=1)


# ---------------------------------------------------------------------------
# read_string round-trips
# ---------------------------------------------------------------------------


def test_read_string_polynomial_only():
    s = "dimension 2\npolynomial\nx0x0: 0.5\nx1x1: 0.5\n"
    polynomial, exponential = read_string(s)
    assert polynomial == {(2, 0): 0.5, (0, 2): 0.5}
    assert exponential == []


def test_read_string_exponential_only():
    s = "dimension 2\nexponential\n1.0, 0.5, [0.0, 0.0]\n"
    polynomial, exponential = read_string(s)
    assert polynomial == {}
    assert len(exponential) == 1
    lincoeff, width, pos = exponential[0]
    assert lincoeff == pytest.approx(1.0)
    assert width == pytest.approx(0.5)
    assert pos == pytest.approx([0.0, 0.0])


def test_read_string_both_sections():
    s = "dimension 2\n" "polynomial\nx0x0: 0.5\n" "exponential\n1.0, 0.5, [0.0, 0.0]\n"
    polynomial, exponential = read_string(s)
    assert (2, 0) in polynomial
    assert len(exponential) == 1


def test_read_string_no_sections_raises():
    with pytest.raises(ValueError, match="No 'polynomial' or 'exponential'"):
        read_string("dimension 2\nx0x0: 0.5\n")


# ---------------------------------------------------------------------------
# obtain_polynomial_squared
# ---------------------------------------------------------------------------


def test_squared_single_monomial():
    # V = x0  →  V² = x0²
    squared = obtain_polynomial_squared({(1,): 1.0})
    assert squared == {(2,): pytest.approx(1.0)}


def test_squared_constant():
    # V = 3  →  V² = 9
    squared = obtain_polynomial_squared({(0,): 3.0})
    assert squared == {(0,): pytest.approx(9.0)}


def test_squared_1d_binomial():
    # V = a·x + b  →  V² = a²·x² + 2ab·x + b²
    a, b = 2.0, 3.0
    squared = obtain_polynomial_squared({(1,): a, (0,): b})
    assert squared.get((2,), 0.0) == pytest.approx(a**2)
    assert squared.get((1,), 0.0) == pytest.approx(2 * a * b)
    assert squared.get((0,), 0.0) == pytest.approx(b**2)


def test_squared_2d_sum():
    # V = x0 + x1  →  V² = x0² + 2·x0·x1 + x1²
    squared = obtain_polynomial_squared({(1, 0): 1.0, (0, 1): 1.0})
    assert squared.get((2, 0), 0.0) == pytest.approx(1.0)
    assert squared.get((0, 2), 0.0) == pytest.approx(1.0)
    assert squared.get((1, 1), 0.0) == pytest.approx(2.0)


def test_squared_result_is_non_negative_for_real_polynomial():
    # For a real polynomial with positive coefficients, all squared coefficients ≥ 0
    polynomial = {(2, 0): 0.5, (0, 2): 0.5}
    squared = obtain_polynomial_squared(polynomial)
    for val in squared.values():
        assert val >= 0.0


# ---------------------------------------------------------------------------
# insert_constants
# ---------------------------------------------------------------------------


def test_insert_constants_substitutes_value_in_poly():
    new_poly, _ = insert_constants(["x0x0: omega"], [], ["omega: 0.5"])
    assert new_poly == ["x0x0: 0.5"]


def test_insert_constants_substitutes_value_in_exp():
    _, new_exp = insert_constants([], ["1.0, alpha, [0.0]"], ["alpha: 0.25"])
    assert new_exp == ["1.0, 0.25, [0.0]"]


def test_read_string_with_constants_section():
    s = "dimension 1\n" "polynomial\nx0x0: omega\n" "constants\nomega: 0.5\n"
    polynomial, _ = read_string(s)
    assert polynomial == {(2,): 0.5}
