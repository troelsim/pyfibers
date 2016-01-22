from pyfibers.atoms import CesiumAtom
from pyfibers.fibers import LeKienFiber
from pyfibers.coupling.base import RadCouplingTensor
import numpy as np
from matplotlib import pyplot as plt

fiber = LeKienFiber(1.45, 1, 1, 1.55)

tensor = RadCouplingTensor(fiber)

Rs = np.linspace(1, 50, 101)

rates = []

nu = 0

for R in Rs:
    print R
    atom = CesiumAtom(R, 0, 0)
    rates.append(tensor.get_coupling(atom, atom, nu, nu))

plt.plot(Rs, rates)
plt.show()
