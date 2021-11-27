import numpy as np
A = np.array([
    [1, 1, 1, 1, 1],
    [-2, -1, 1, 2, 3],
    [2, 1 / 2, 1 / 2, 2, 9 / 2],
    [-4 / 3, -1 / 6, 1 / 6, 4 / 3, 9 / 2],
    [2 / 3, 1 / 24, 1 / 24, 2 / 3, 27 / 8]
], dtype=np.float64)
b = np.array([[0, 0, 0, 1, 0]]).T

print(np.linalg.solve(A, b).T)