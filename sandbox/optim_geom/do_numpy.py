import numpy as np
from numpy.linalg import eigh

Natoms = 2
H = np.eye(Natoms*3)*70
H[1,1] = 22.0
H[2,3] = 11.0
H[3,2] = 11.0
print(H)

f = np.array([ [3.0, 2.0, 3.0],
               [3.1, 2.1, 4.0] ])
print(f)

f = f.reshape(-1)
omega, V = eigh(H)
print(omega)
dr = np.dot( V, np.dot(f, V) / np.fabs(omega) ).reshape((-1, 3))
print(dr)

steplengths = (dr**2).sum(1)**0.5
print(steplengths)
