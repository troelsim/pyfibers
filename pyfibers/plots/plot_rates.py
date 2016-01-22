from pyfibers.atoms import CesiumAtom
from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import  LeKienRadMode
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt


fiber = LeKienFiber(1.45, 1, 1, 1.55)

points = 101
Vs = np.sqrt(1.45**2-1)*np.linspace(0.1, 12.2, points)
V=1.55

Rs = np.linspace(1.0, 5, points)
m_max = 7
couplings = np.zeros((points, 2*m_max+1))
errs = np.zeros((points, 2*m_max+1))

atom = CesiumAtom(1,0,0)
Me = 2

#for V in Vs:
#    print V
for i, R in enumerate(Rs):
    fiber = LeKienFiber(1.45, 1, 1, V)
    pos = (R,0,0)
    pol = 1
    f = 1
    rates = []
    for m in range(-m_max, m_max+1):
        def integrand(rb):
            U = np.sqrt((fiber.n*fiber.rk)**2 - rb**2)
            mode = LeKienRadMode(fiber, m, U, pol, f)
            norm = mode.norm()
            rate = 0
            vec = np.array([
                mode.e_r(*pos),
                mode.e_phi(*pos),
                mode.e_z(*pos),
            ])
            for dv in [atom.dipole_vector(Me, i) for i in [Me-1, Me, Me+1]]:
                rate += abs(np.dot(vec, dv))**2
            return rate/norm
        res, err = quad(integrand, 0.01, fiber.nc*fiber.rk-0.01)
        rates.append(res)
        couplings[i, m] = res
        errs[i, m] = err

    print rates
    print "Max: %e, min: %e, mean: %e, median: %e" % (max(rates), min(rates), np.average(rates), np.median(rates))
    avg_over_median = np.average(rates)/np.median(rates)
    if avg_over_median < 0.2 or avg_over_median > 5:
        print "AVERAGE AND MEDIAN DIFFER BY MORE THAN ORDER OF MAGNITUDE"

for i in range(-m_max, m_max+1):
    plt.plot(Rs, np.sum(couplings[:,:i+1], 1), label=i)
plt.legend()
plt.show()
