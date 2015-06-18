from numpy import pi, sqrt
from numpy.linalg import norm
import numpy as np
from pyfibers.atoms import CesiumAtom
from pyfibers.fibers import LeKienFiber
from pyfibers.coupling.fortran import FortranCouplingTensor
from pyfibers.atoms import SimpleAtom
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt


rho = 200.e-9
wavelength = 854e-9
k = 2*pi/wavelength

#k = k*rho
#rho = 1.0

n = 1.45
nc = 1.0

fiber = LeKienFiber(n, nc, rho, rho*k*sqrt(n**2-nc**2))

Fe = CesiumAtom.Fe
Fg = CesiumAtom.Fg

Me = 3
#@np.vectorize
def rate(R):
    atom = CesiumAtom(R, 0, 0)
    Mgs = set(range(Me-1,Me+2)) & set(range(-Fg, Fg+1))
    sum = 0
    dipole_strength_total = 0

    ft = FortranCouplingTensor(fiber)
    coupling_matrix = ft.get_coupling_matrix(atom, atom)

    for Mg in Mgs:
        dipole_vector = atom.dipole_vector(Me, Mg)
        dipole_strength = norm(dipole_vector)**2

        coupling = (dipole_vector.transpose()*coupling_matrix*dipole_vector.conjugate()).item(0,0)
        dipole_strength_total += dipole_strength
        sum += coupling

    return sum/dipole_strength_total



from multiprocessing import Pool
p=Pool(5)
Rs = np.linspace(1.0,10,51)
es = p.map(rate, Rs)

plt.plot(Rs, es)
plt.semilogy()
plt.show()
print "Saving..."
#plt.savefig('/Users/troelsim/ftest.png')
print "Done!"
print rate(3, 1)
