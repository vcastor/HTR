import numpy as np
from numpy import linalg as LA
from scipy.linalg import eigh
#import sympy as sp

F = np.array([[-2.61253, -1.43358], [-1.43358, -1.73439]])
S = np.array([[1., 0.53707], [0.53707, 1.]])

E,C = eigh(F, S, eigvals_only=False)
#C = eigh(S, eigvals_only=False)

print(E)
print(C)
#E,C = eigh(F, S, eigvals_only=False)
