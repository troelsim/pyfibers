from numpy.linalg import norm
from scipy.integrate import quad
from pyfibers.atoms import CesiumAtom
from pyfibers.constants import SOL
from pyfibers.fibers import LeKienFiber
from numpy import sqrt, pi, vectorize
from pyfibers.modes.lekien import LeKienRadMode
import numpy as np
from matplotlib import pyplot as plt

atom = CesiumAtom(1,0,0)  # Cesium atom at the fiber surface, as in LK

rho = 200.e-9
wavelength = 854e-9
k = 2*pi/wavelength
n = 1.45
nc = 1.0

fiber = LeKienFiber(n, nc, rho, rho*k*sqrt(n**2-nc**2))

Fe = atom.Fe
Fg = atom.Fg

@np.vectorize
def rate(Me, R):
    atom = CesiumAtom(R,0,0)
    # There is a nonzero coupling for |Me-Mg| <= 1
    Mgs = set(range(Me-1,Me+2)) & set(range(-Fg, Fg+1))
    sum = 0
    dipole_strength_total = 0
    for Mg in Mgs:
        dipole_vector = atom.dipole_vector(Me, Mg)
        dipole_strength = norm(dipole_vector)**2

        def integrand(beta):
            U = sqrt(fiber.rk**2*fiber.n**2 - fiber.rho**2*beta**2)
            m = LeKienRadMode(fiber, 1, U)
            result = 0
            field_norm = m.norm()
            for nu in range(-5, 6):
                m.nu = nu
                e_vector = np.matrix(m.e_vector(atom)).transpose()
                coupling = norm(dipole_vector.transpose()*e_vector)**2
                result += coupling
            return result/field_norm

        delta = 1e-3
        res, eps = quad(integrand, delta, k*(1-delta))
        print res, eps
        sum += res*6*pi*SOL/k**2
        dipole_strength_total += dipole_strength
    return sum

print rate(3,1)

@np.vectorize
def nn(beta):
    U = sqrt(fiber.rk**2*fiber.n**2 - fiber.rho**2*beta**2)
    m = LeKienRadMode(fiber, 1, U)
    return m.norm()

points = np.linspace(1,2,31)
#plt.plot(points, nn(points))
plt.plot(points, rate(3,points))
plt.show()