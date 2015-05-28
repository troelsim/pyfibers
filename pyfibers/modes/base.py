from numpy import array
from .util import vectorize_method

class FiberMode(object):
    """
    Mode base class, for guided as well as radiative modes
    """
    @vectorize_method
    def e_r(self, R, phi, z):
        """
        Delegates the calculation to e_r_core and e_r_cladding
        :param R:
        :param phi:
        :param z:
        :return:
        """
        if R <= 1:
            result = self.e_r_core(R, phi, z)
        else:
            result = self.e_r_cladding(R, phi, z)
        return result*self.implicit(R, phi, z)

    @vectorize_method
    def e_phi(self, R, phi, z):
        """
        Delegates the calculation to e_phi_core and e_phi_cladding
        :param R:
        :param phi:
        :param z:
        :return:
        """
        if R <= 1:
            result = self.e_phi_core(R, phi, z)
        else:
            result = self.e_phi_cladding(R, phi, z)
        return result*self.implicit(R, phi, z)

    @vectorize_method
    def e_z(self, R, phi, z):
        """
        Delegates the calculation to e_z_core and e_z_cladding
        :param R:
        :param phi:
        :param z:
        :return:
        """
        if R <= 1:
            result = self.e_z_core(R, phi, z)
        else:
            result = self.e_z_cladding(R, phi, z)
        return result*self.implicit(R, phi, z)

    @vectorize_method
    def e_vector(self, atom):
        """
        Returns the field vector in (R, phi, z) basis at an atom's position
        :param atom: Atom
        :return: NumPy array (column)
        """
        return array([
            self.e_r(atom.R, atom.phi, atom.z),
            self.e_phi(atom.R, atom.phi, atom.z),
            self.e_z(atom.R, atom.phi, atom.z)
        ]).transpose()

    def e_r_core(self, R, phi, z):
        pass

    def e_r_cladding(self, R, phi, z):
        pass

    def e_phi_core(self, R, phi, z):
        pass

    def e_phi_cladding(self, R, phi, z):
        pass

    def e_z_core(self, R, phi, z):
        pass

    def e_z_cladding(self, R, phi, z):
        pass

    def implicit(self, R, phi, z):
        return 1.0 + 0.0j

    def discrete_modes(self, *args, **kwargs):
        pass
