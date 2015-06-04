"""
NOTE: Not a unit test
"""
from numpy.linalg import norm

from pyfibers.modes.lekien import LeKienGuidedMode
from pyfibers.fibers import LeKienFiber
from pyfibers.atoms import CesiumAtom
from pyfibers.constants import *
from numpy import pi, sqrt, vectorize
import numpy as np

from matplotlib import pyplot as plt

atom = CesiumAtom(1,0,0)  # Cesium atom at the fiber surface, as in LK

rho = 200.e-9
wavelength = 854e-9
k = 2*pi/wavelength
n = 1.45
nc = 1.0

fiber = LeKienFiber(n, nc, rho, rho*k*sqrt(n**2-nc**2))

mode = LeKienGuidedMode(fiber, 1, 1)  # Only the +,+ mode is necessary because of symmetry


Fe = atom.Fe
Fg = atom.Fg

@vectorize
def decay_rate(e, R):
    atom = CesiumAtom(R, 0, 0)
    e_vector = mode.e_vector(atom)
    d_vector = sum([atom.dipole_vector(e, g) for g in range(-4,5)])
    coupling = 0
    dipole_strength = 0
    for g in range(-Fg,Fg+1):
        d_vector = atom.dipole_vector(e, g)
        coupling += (np.abs(e_vector.transpose() * d_vector)**2).item(0)
        dipole_strength += norm(d_vector)**2

    return mode.betaprime*SOL*coupling/mode.closed_norm/dipole_strength * 6*pi/k**2

Rs = np.linspace(1,3,101)

for e in range(-Fe,Fe+1):
    rates = decay_rate(e, Rs)
    plt.plot(Rs, rates)

plt.show()