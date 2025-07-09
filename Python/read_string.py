import re
import sys
example_string="""
dimension 4
polynomial
x0x1: 2
x2x3: 1
x0x2: 1
x1x3: -1
"""

def read_string(string: str):
    """
    Read a string containing polynomial coefficients and return a dictionary.
    
    The string should contain lines with the format:
    "x0: 1", "x1x2: 2", etc.
    
    Returns:
        dict: A dictionary with keys as variable names and values as coefficients.
    """
    result={}
    lines=string.splitlines()
    while lines[0].strip() == "": # Skip empty lines
        lines.pop(0)
    dimension= int(lines[0].strip().split()[1])
    typex= lines[1].strip()
    poly_list= []
    print("Dimension:", dimension
          , "Type:", typex)
    if typex != "polynomial":
        raise ValueError(f"Expected 'polynomial' type, got '{typex}'")
    else:
        for line_idx, raw in enumerate(lines[2:]):       # outer index ≠ 'i'
            poly, val = raw.strip().split(":")
            line_vals = []
            for d in range(dimension):                   # use 'd' not 'i'
                hits = re.findall(fr"x{d}(?!\d)", poly)  # exact matches only
                line_vals.append(len(hits))
            if line_vals in poly_list:
                raise ValueError(f"Duplicate polynomial found at line {line_idx + 3}: {raw.strip()}")
            poly_list.append(line_vals)
            result[tuple(line_vals)] = float(val.strip())
    return result
polynomial=read_string(example_string)
def obtain_polynomial_squared(polynomial):
    """
    Obtain the squared polynomial coefficients from the given polynomial.
    
    Args:
        polynomial (list): List of lists containing variable counts and coefficients.
        
    Returns:
        list: List of lists with squared polynomial coefficients.
    """
    result={}
    iterator=list(polynomial.items())
    for i, (vars1, coeff1) in enumerate(iterator):
        for j, (vars2, coeff2) in enumerate(iterator):
            if i<j:
                continue
            duplicator= 1 if i==j else 2 # Take into account double-counting of non-diagonal terms
            new_vars = [v1 + v2 for v1, v2 in zip(vars1, vars2)]
            new_coeff = coeff1 * coeff2* duplicator
            result[tuple(new_vars)] = result.get(tuple(new_vars), 0) + new_coeff
    return result
polynomial_squared = obtain_polynomial_squared(polynomial)
print("Original Polynomial Coefficients:")
print(polynomial)
print("Squared Polynomial Coefficients:")
print(polynomial_squared)