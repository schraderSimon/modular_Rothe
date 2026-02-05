#!/usr/bin/env python3

import sys

sys.path.append(".")

from libraries.file_handling import load_rothe_file

# Load the data and see what keys it has
data = load_rothe_file("Hydrogen_50_21", "none", path="./wave_function_data")

print("Keys in the data:")
for key in data.keys():
    if isinstance(data[key], list):
        print(f"  {key}: list with {len(data[key])} elements")
        if len(data[key]) > 0:
            print(f"    First element type: {type(data[key][0])}")
            if hasattr(data[key][0], "shape"):
                print(f"    First element shape: {data[key][0].shape}")
    else:
        print(f"  {key}: {type(data[key])}")
        if hasattr(data[key], "shape"):
            print(f"    Shape: {data[key].shape}")

print("\nFirst few time points:")
print(data["t"][:5] if len(data["t"]) > 5 else data["t"])
