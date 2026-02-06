example_string_6D = """
    dimension 6
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x2x2: 0.5
    x3x3: 0.5
    x4x4: 0.5
    x5x5: 0.5
    x0x0x1: 0.111803
    x1x1x2: 0.111803
    x2x3x3: 0.111803
    x3x4x4: 0.111803
    x4x5x5: 0.111803
    x1x1x1: -0.03726766666
    x2x2x2: -0.03726766666
    x3x3x3: -0.03726766666
    x4x4x4: -0.03726766666
    x5x5x5: -0.03726766666
    """

example_string_4D = """
    dimension 4
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x2x2: 0.5
    x3x3: 0.5
    x0x0x1: 0.111803
    x1x1x2: 0.111803
    x2x3x3: 0.111803
    x1x1x1: -0.03726766666
    x2x2x2: -0.03726766666
    x3x3x3: -0.03726766666
    """
example_string_2D = """
    dimension 2
    polynomial
    x0x0: 0.5
    x1x1: 0.5
    x0x0x1: 0.111803
    x1x1x1: -0.03726766666
    x0x0x0x0: 0.0007812444255625
    x1x1x1x1: 0.0007812444255625
    x0x0x1x1: 0.001562488851125
    """
example_string_1D = """
    dimension 1
    polynomial
    x0x0: 0.5
    x0x0x0x0: 0.0007812444255625
    """
henonheiles_strings = {}
henonheiles_strings["1D"] = example_string_1D
henonheiles_strings["2D"] = example_string_2D
henonheiles_strings["4D"] = example_string_4D
henonheiles_strings["6D"] = example_string_6D
