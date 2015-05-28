from math import sqrt
from sympy.physics.wigner import wigner_3j, wigner_6j

class Atom(object):
    def __init__(self, R, phi, z):
        self.R = R
        self.phi = phi
        self.z = z


class CesiumAtom(Atom):
    reduced_dipole_element = 1.0

    # Upper (F=5) quantum numbers
    Le = 1
    Se = 0.5
    Ie = 3.5
    Je = 1.5
    Fe = 5

    # Lower (F=4) quantum numbers
    Lg = 0
    Sg = 0.5
    Ig = 3.5
    Jg = 0.5
    Fg = 4

    def dipole(self, Me, Mg, q):
        """
        Find the q-spherical component of the dipole matrix element between Me and Mg
        :param M1:
        :param M2:
        :param q:
        :return:
        """
        return (
            (-1)**(self.Ig + self.Je - Me) * self.reduced_dipole_element *
            sqrt((2*self.Fg+1)*(2*self.Fe+1)) *
            float(wigner_6j(self.Je, self.Fe, self.Ig, self.Fg, self.Jg, 1)) *
            float(wigner_3j(self.Fg, 1, self.Fe, Mg, q, -Me))
        )


class SimpleAtom(CesiumAtom):
    """
    An atom consisting of only the leftmost levels of the Cesium Atom
    """

    def dipole(self, Me, q):
        """
        The dipole moment for Me = -1, 0, 1
        :param Me:
        :param q:
        :return:
        """
        assert(Me in [-1, 0, 1])
        Mg = 4
        return super(SimpleAtom, self).dipole(Me+4, Mg, q)