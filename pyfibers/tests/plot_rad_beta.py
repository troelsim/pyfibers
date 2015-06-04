# Plot radial fields as a function of beta

from matplotlib import pyplot as plt
from numpy.linalg import norm
from pyfibers.atoms import CesiumAtom
from pyfibers.fibers import LeKienFiber
from numpy import pi, sqrt
from pyfibers.modes.lekien import LeKienRadMode
import numpy as np

atom = CesiumAtom(1.1,0,0)  # Cesium atom at the fiber surface, as in LK

rho = 200.e-9
wavelength = 854e-9
k = 2*pi/wavelength
n = 1.45
nc = 1.0

fiber = LeKienFiber(n, nc, rho, rho*k*sqrt(n**2-nc**2))

@np.vectorize
def e_norm(beta):
    U = sqrt(fiber.rk**2*fiber.n**2 - fiber.rho**2*beta**2)
    m = LeKienRadMode(fiber, 1, U)
    return norm(m.e_vector(atom))/sqrt(m.norm())

betas = np.linspace(0, k, 101)
plt.plot(betas, e_norm(betas))
plt.show()