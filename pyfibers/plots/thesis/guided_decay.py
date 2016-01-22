from pyfibers.shortcuts import lkfiber as fiber
from pyfibers.coupling.base import GuidedCouplingTensor
from pyfibers.atoms import SimpleAtom
import numpy as np
from matplotlib import pyplot as plt

ct = GuidedCouplingTensor(fiber)

# Two atoms on line, variable distance

points = 201
zs = np.linspace(0, 10, points)

atom1 = SimpleAtom(1, 0, 0)


rates = np.zeros((6, points), dtype='complex')

for point, z in enumerate(zs):
    master_matrix = np.zeros((6, 6), dtype='complex')
    for i in range(6):
        for j in range(6):
            # Translate matrix indices to quantu numbers
            nui = 1 if i % 3 == 0 else (-1 if i % 3==1 else 0)
            nuj = 1 if j % 3 == 0 else (-1 if j % 3==1 else 0)
            atom2 = SimpleAtom(1, 0, z)
            ct.get_coupling_matrix(atom1, atom2)
            master_matrix[i, j] = ct.get_coupling(atom1, atom2, nui, nuj)
    print "Is it Hermitian?"
    print master_matrix - master_matrix.transpose().conjugate()
    values, vectors = np.linalg.eig(master_matrix)
    print "Eigenvalues: ", values
    rates[:, point] = values

for line in range(6):
    plt.plot(zs, rates[line,:])

plt.show()