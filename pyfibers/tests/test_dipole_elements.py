from pyfibers.atoms import CesiumAtom
from numpy.linalg import norm
import numpy as np
import sys

atom = CesiumAtom(1,0,0)
Fe = atom.Fe
Fg = atom.Fg

cp = np.zeros((2*Fe+1, 2*Fg+1))

for Me in range(-Fe, Fe+1):
    for Mg in range(-Fg, Fg+1):
        n = norm(atom.dipole_vector(Me, Mg))
        cp[Me+Fe, Mg+Fg] = n
        sys.stdout.write("%f\t" % n)
    print ""


for i in range(2*Fg+1):
    print sum(cp[:,i])