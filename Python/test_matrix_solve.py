import sys
import numpy as np
from scipy.linalg import solve
import time

n = int(sys.argv[1])
A = np.random.randn(n, n) + 1j * np.random.randn(n, n)
A = A @ A.conj().T  # Makes it Hermitian positive-definite
b = np.random.randn(n) + 1j * np.random.randn(n)

start = time.time()
x = solve(A, b, assume_a='pos')
end = time.time()

print(f"t({n}) = {end - start:.6f} seconds")
