import numpy as np
from scipy.integrate import quad
from pyfibers.modes.lekien import LeKienRadMode
from pyfibers.constants import HBAR


class CouplingTensor(object):
    def spherical_basis(self, phi):
        return np.matrix([
            [-np.exp(-1j*phi), 1j*np.exp(1j*phi), 0         ],
            [np.exp(-1j*phi),  1j*np.exp(1j*phi), 0         ],
            [0,                0,                 np.sqrt(2)]
        ], dtype='c16')/np.sqrt(2)

    def get_coupling(self, atom_j, atom_k, nu_j, nu_k):
        pass

    def get_coupling_matrix(self, atom_j, atom_k):
        pass


class RadCouplingTensor(CouplingTensor):

    def __init__(self, fiber):
        self.fiber = fiber
        self.mode_class = fiber.rad_mode_class
        self.calls = 0

    def get_coupling(self, atom_j, atom_k, nu_j, nu_k):
        Uj = self.spherical_basis(atom_j.phi)
        Uk = self.spherical_basis(atom_k.phi)

        gamma_0 = atom_j.emission_rate(nu_j, 0, self.fiber.rk/self.fiber.rho)

        def integrand(rb):
            U = self.mode_class.get_U(self.fiber, rb)
            coupling = 0.+0j
            for mode in self.mode_class.discrete_modes(self.fiber, U=U, m_max=4):
                e_j = np.matrix(mode.e_vector(atom_j)).transpose()
                e_k = np.matrix(mode.e_vector(atom_k)).transpose()

                d_j = np.matrix([
                    atom_j.dipole(nu_j, 0, 1),
                    atom_j.dipole(nu_j, 0, -1),
                    atom_j.dipole(nu_j, 0, 0),
                ]).transpose()

                d_k = np.matrix([
                    atom_k.dipole(nu_k, 0, 1),
                    atom_k.dipole(nu_k, 0, -1),
                    atom_k.dipole(nu_k, 0, 0),
                ]).transpose()

                tensor_matrix = e_j*e_k.transpose().conjugate()

                coupling += - d_j.transpose() * tensor_matrix * d_k.conjugate() / mode.norm() / gamma_0* 1/(2*np.pi*HBAR)
            self.calls += 1
            return coupling[0, 0]

        res, err = quad(lambda x: np.real(integrand(x)), 0.01*self.fiber.rk, self.fiber.nc*self.fiber.rk*0.99)
        try:
            print "Relative error in integral: %e" % err/res
        except TypeError:
            print res, err
        return res

    def get_coupling_matrix(self, atom_j, atom_k):
        """
        The coupling between j and k upper levels, in matrix form
        :param atom_j: One atom
        :param atom_k: Another atom
        :return:
        """
        matrix = np.matrix(np.zeros((3,3)))

        def integrand(beta, Mj, Ml):
            U = np.sqrt(self.fiber.rk**2*self.fiber.n**2 - self.fiber.rho**2*beta**2)
            coupling = 0.+0j
            for mode in self.mode_class.discrete_modes(self.fiber, U=U):
                e_j = np.matrix(mode.e_vector(atom_j)).transpose()
                e_k = np.matrix(mode.e_vector(atom_k)).transpose()

                tensor_matrix = e_j*e_k.transpose().conjugate()
                dj = np.matrix(atom_j.dipole_vector(Mj, 0))
                dk = np.matrix(atom_k.dipole_vector(Ml, 0))
                coupling += (dj.transpose() * tensor_matrix * dk.conjugate() / mode.norm()).item(0)
            return np.real(coupling)

        for Mj in [1, -1, 0]:
            for Ml in [1, -1, 0]:
                res, tol = quad(integrand, 0, self.fiber.rk/self.fiber.rho, args=(Mj, Ml))
                print "%d, %d: %e (+/- %e)" % (Mj, Ml, res, tol)
                matrix[Mj+1,Ml+1] += res

        return matrix


class GuidedCouplingTensor(CouplingTensor):
    def __init__(self, fiber):
        self.fiber = fiber
        self.mode_class = fiber.guided_mode_class

    def get_coupling(self, atom_j, atom_k, nu_j, nu_k):
        Uj = self.spherical_basis(atom_j.phi)
        Uk = self.spherical_basis(atom_k.phi)

        gamma_0 = atom_j.emission_rate(nu_j, 0, self.fiber.rk/self.fiber.rho)
        coupling = 0.+0j
        for mode in self.mode_class.discrete_modes(self.fiber):
            e_j = np.matrix(mode.e_vector(atom_j)).transpose()
            e_k = np.matrix(mode.e_vector(atom_k)).transpose()

            d_j = np.matrix([
                atom_j.dipole(nu_j, 0, 1),
                atom_j.dipole(nu_j, 0, -1),
                atom_j.dipole(nu_j, 0, 0),
            ]).transpose()

            d_k = np.matrix([
                atom_k.dipole(nu_k, 0, 1),
                atom_k.dipole(nu_k, 0, -1),
                atom_k.dipole(nu_k, 0, 0),
            ]).transpose()

            tensor_matrix = e_j*e_k.transpose().conjugate()

            coupling += d_j.transpose() * tensor_matrix * d_k.conjugate() / mode.norm() / gamma_0 * 1/(2*np.pi*HBAR)

        return coupling[0, 0]

    def get_coupling_matrix(self, atom_j, atom_k):
        """
        The coupling between j and k upper levels, in matrix form
        :param atom_j: One atom
        :param atom_k: Another atom
        :return:
        """
        matrix = np.matrix(np.zeros((3,3)))
        for mode in self.mode_class.discrete_modes(self.fiber):
            e_j = np.matrix(mode.e_vector(atom_j)).transpose()
            e_k = np.matrix(mode.e_vector(atom_k)).transpose()

            tensor_matrix = e_j*e_k.transpose().conjugate()

            for Mj in [1, -1, 0]:
                for Ml in [1, -1, 0]:
                    dj = np.matrix(atom_j.dipole_vector(Mj, 0))
                    dk = np.matrix(atom_k.dipole_vector(Ml, 0))
                    coupling = (dj.transpose() * tensor_matrix * dk.conjugate() / mode.norm()).item(0)
                    matrix[Mj+1,Ml+1] += coupling

        return matrix

if __name__ == '__main__':
    from pyfibers.atoms import SimpleAtom
    from pyfibers.fibers import LeKienFiber
    f = LeKienFiber(1.45,1,1,1)
    a1 = SimpleAtom(1,0,0)
    a2 = SimpleAtom(1,0,1)
    t = GuidedCouplingTensor(f)
    tr = RadCouplingTensor(f)
    print t.get_coupling_matrix(a1, a2)
    print tr.get_coupling_matrix(a1, a2)


def get_coupling_matrix(fiber, r, zs, phis):
    # Get the (guided) coupling matrix for a given fiber
    assert len(zs) == len(phis)
    num_atoms = len(zs)
    from pyfibers.modes.lekien import LeKienGuidedMode
    mode = LeKienGuidedMode(fiber, pol=1, f=1)
    e_r = -1j*mode.e_r(r, 0, 0)
    e_phi = mode.e_phi(r, 0, 0)
    e_z = mode.e_z(r, 0, 0)

    def spherical_basis(phi):
        return np.array([
            [-np.exp(1j*phi), np.exp(1j*phi), 0],
            [-1j*np.exp(-1j*phi), -1j*np.exp(-1j*phi), 0],
            [0, 0, np.sqrt(2)]
        ])/np.sqrt(2)

    coupling_matrix = np.zeros((3*num_atoms, 3*num_atoms), dtype='complex')
    def sub_matrix(l, j):
        dz = zs[l]-zs[j]
        phil = phis[l]
        phij = phis[j]
        dphi = phis[l]-phis[j]
        b = mode.rb/fiber.rho
        cosb = np.cos(b*dz)
        cosphi = np.cos(dphi)
        sinb = np.sin(b*dz)
        sinphi = np.sin(dphi)
        mat = np.zeros((3,3), dtype='complex')
        return np.dot(spherical_basis(phil), np.dot(np.array([
            [e_r**2*cosb*cosphi, -e_r*e_phi*cosb*sinphi, -e_r*e_z*sinb*cosphi],
            [e_phi*e_r*cosb*sinphi, e_phi**2*cosb*cosphi, e_phi*e_z*sinb*sinphi],
            [e_z*e_r*sinb*cosphi, e_z*e_phi*sinb*sinphi, e_z**2*cosb*cosphi]
        ]), spherical_basis(phij).transpose().conjugate()))

    for l in range(num_atoms):
        for j in range(num_atoms):
            coupling_matrix[3*l:3*(l+1), 3*j:3*(j+1)] = sub_matrix(l, j)

    return coupling_matrix


def get_e_matrix(fiber, r, zs, phis, Vmin=0.1, Vmax=2.405):
    V0 = fiber.V
    N = len(zs)
    from pyfibers.fibers import LeKienFiber
    def matrix_el(i, j):
        def integrand(V):
            fib = LeKienFiber(fiber.n, fiber.nc, fiber.rho, V)
            return get_coupling_matrix(fib, r, zs, phis)[i,j]*V/(V-V0)
        def fold(x):
            # fold the functino around V0
            return integrand(V0+x) + integrand(V0-x)
        eps=0.1
        i1, p1 = quad(integrand, Vmin, V0-eps)
        print i1, p1
        i2, p2 = quad(fold, 0, eps)
        print i2, p2
        i3, p3 = quad(integrand, V0+eps, Vmax)
        print i3, p3
        print i1 + i2 + i3
        return i1 + i2 + i3
    res = np.zeros((3*N, 3*N), dtype='complex')
    for i in range(0, 3*N):
        for j in range(0, 3*N):
            res[i,j] = matrix_el(i, j)
    return res

def get_e_eigenvalue(fiber, r, dz):
    mat = get_e_matrix(fiber, r, [0, dz], [0, 0], Vmin=0.76, Vmax=2.405)
    eigs, eigvecs = np.linalg.eigh(mat)
    return eigs

