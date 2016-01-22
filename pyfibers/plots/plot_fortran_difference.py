
import numpy as np
from matplotlib import pyplot as plt
from ffibers.fibers import fibers

from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import LeKienRadMode

V = 1.55
n = 1.45
nc = 1.
rk = V/np.sqrt(n**2-nc**2)

rb = 0.5*rk

nu = 1
U = np.sqrt(n**2*rk**2-rb**2)

pol = 1
f = 1

Rs = np.linspace(0.1, 10, 101)

fibers.init_fiber(n, nc, 1., V)
fibers.init_mode(nu, U, pol, f)

fiber = LeKienFiber(n, nc, 1, V)
mode = LeKienRadMode(fiber, nu, U, pol, f)

def discgen(mmax):
    for m in range(-mmax, mmax+1):
        for pol_ in [-1, 1]:
            for f_ in [-1, 1]:
                yield (m, pol_, f_)

for nu, pol, f in discgen(2):
    for dim in ['r', 'phi', 'z']:
        f = getattr(fibers, "e_%s" % dim)
        esf = np.vectorize(f)(Rs, 0, 0)/fibers.norm()
        print esf
        p = getattr(mode, "e_%s" % dim)
        esp = p(Rs, 0, 0)/mode.norm()
        plt.plot(Rs, esf, 'k', Rs, esp, ':k')
        plt.legend()
    plt.show()
