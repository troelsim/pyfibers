import numpy as np
from scipy.integrate import quad
from pyfibers.modes.lekien import LeKienRadMode


class CouplingTensor(object):
    def spherical_basis(self, phi):
        return np.matrix([
            [-np.exp(-1j*phi), 1j*np.exp(1j*phi), 0         ],
            [np.exp(-1j*phi),  1j*np.exp(1j*phi), 0         ],
            [0,                0,                 np.sqrt(2)]
        ], dtype='c16')/np.sqrt(2)


class RadCouplingTensor(object):

    def __init__(self, fiber):
        self.fiber = fiber
        self.mode_class = fiber.rad_mode_class
        self.calls = 0

    def get_coupling(self, atom_j, atom_k, nu_j, nu_k):
        Uj = self.spherical_basis(atom_j.phi)
        Uk = self.spherical_basis(atom_k.phi)

        def integrand(U):
            coupling = 0.+0j
            for mode in self.mode_class.discrete_modes(self.fiber, U=U, m_max=4):
                e_j = np.matrix(mode.e_vector(atom_j)).transpose()
                e_k = np.matrix(mode.e_vector(atom_k)).transpose()

                d_j = np.matrix([
                    atom_j.dipole(nu_j, 1),
                    atom_j.dipole(nu_j, -1),
                    atom_j.dipole(nu_j, 0),
                ]).transpose()

                d_k = np.matrix([
                    atom_k.dipole(nu_k, 1),
                    atom_k.dipole(nu_k, -1),
                    atom_k.dipole(nu_k, 0),
                ]).transpose()

                tensor_matrix = e_j*e_k.transpose().conjugate()

                coupling += d_j.transpose() * tensor_matrix * d_k.conjugate() / mode.norm()
            self.calls += 1
            return coupling[0, 0]

        return quad(lambda x: np.real(integrand(x)), self.mode_class.U_min(self.fiber), self.mode_class.U_max(self.fiber))


class GuidedCouplingTensor(object):
    def __init__(self, fiber):
        self.fiber = fiber
        self.mode_class = fiber.rad_mode_class
