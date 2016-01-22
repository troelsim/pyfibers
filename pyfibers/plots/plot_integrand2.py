"""
Plotting the integrand, without the Fortran module
"""


from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import  LeKienRadMode
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt


fiber = LeKienFiber(1.45, 1, 1, 1.55)

points = 11
Vs = np.sqrt(1.45**2-1)*np.linspace(0.1, 12.2, points)
V=1.55

Rs = np.linspace(1.0, 5, points)
m_max = 7
couplings = np.zeros((points, 2*m_max+1))
errs = np.zeros((points, 2*m_max+1))

#for V in Vs:
#    print V

R = 1
fiber = LeKienFiber(1.45, 1, 1, V)
pol = 1
f = 1

rates = []

rbs = np.linspace(0.01, fiber.rk-0.01)
for R in np.linspace(1,5,11):
    print R
    pos = (R,0,0)
    for m in range(-m_max,m_max+1):
        def integrand(rb):
            U = np.sqrt((fiber.n*fiber.rk)**2 - rb**2)
            mode = LeKienRadMode(fiber, m, U, pol, f)
            norm = mode.norm()
            if not norm > 0:
                print "NORM IS %e" % norm
                mode.norm()
            return (
                abs(mode.e_z(*pos))**2 +
                abs(mode.e_r(*pos))**2 +
                abs(mode.e_phi(*pos))**2
            )/norm
        try:
            plt.plot(rbs, np.vectorize(integrand)(rbs), label=m)
        except ZeroDivisionError:
            print "Zero division at m=%d" % m
            pass
    plt.semilogy()
    plt.show()
    plt.legend()

