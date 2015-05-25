import numpy as np
from pyfibers.modes.lekien import LeKienRadMode


class RadCouplingTensor(object):
    def __init__(self, fiber):
        self.fiber = fiber
        self.mode_class = fiber.rad_mode_class

    def spherical_basis(self, phi):
        return np.matrix(
            [-np.exp(-1j*phi), 1j*np.exp(1j*phi), 0         ],
            [np.exp(-1j*phi),  1j*np.exp(1j*phi), 0         ],
            [0,                0,                 np.sqrt(2)]
        )/np.sqrt(2)

    def get_matrix(self, R1, R2, phi1, phi2, z1, z2):
        def integrand(U):
            for mode in self.mode_class.discrete_modes(U):
                e_j = np.array([
                    mode.e_r(R1, phi1, z1),
                    mode.e_phi(R1, phi1, z1),
                    mode.e_z(R1, phi1, z1)
                ]).transpose()  # column
                e_k = np.array([
                    mode.e_r(R2, phi2, z2),
                    mode.e_phi(R2, phi2, z2),
                    mode.e_z(R2, phi2, z2)
                ])  # row

                tensor_matrix = np.dot(e_j, e_k.conjugate())




#    def get_element(self):
