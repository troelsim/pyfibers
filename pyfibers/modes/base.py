from numpy import vectorize

class FiberMode(object):
    """
    Mode base class, for guided as well as radiative modes
    """
    @vectorize
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

    @vectorize
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

    @vectorize
    def e_z(self, R, phi, z):
        """
        Delegates the calculation to e_z_core and e_z_cladding
        :param R:
        :param phi:
        :param z:
        :return:
        """
        if R <= 1:
            return self.e_z_core(R, phi, z)
        else:
            return self.e_z_cladding(R, phi, z)
        return result*self.implicit(R, phi, z)

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
