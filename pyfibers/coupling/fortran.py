from functools import wraps
from .base import CouplingTensor
from ffibers.fibers import fibers
import numpy as np
from pyfibers.constants import *

def init_first(func):
    @wraps(func)
    def inner(obj, *args, **kwargs):
        obj.init_fiber()
        return func(obj, *args, **kwargs)
    return inner


class FortranCouplingTensor(CouplingTensor):
    fiber = None

    def __init__(self, fiber):
        self.fiber = fiber

    def init_fiber(self):
        fiber = self.fiber
        fibers.init_fiber(fiber.n, fiber.nc, fiber.rho, fiber.V)

    @init_first
    def get_coupling_matrix(self, atom_j, atom_k):
        """
        Coupling matrix using external fortran library
        :param atom_j:
        :param atom_k:
        :return:
        """
        coupling_matrix = np.zeros((3, 3), dtype='c16')
        for f in [-1,1]:
            for pol in [-1, 1]:
                for m in [-5, 6]:
                    fibers.init_mode(m, self.fiber.V, pol, f)
                    for i in range(3):
                        for j in range(3):
                            coupling_matrix[i, j] += fibers.coupling(i, j,
                                                                     atom_j.R, atom_j.phi, atom_j.z,
                                                                     atom_k.R, atom_k.phi, atom_k.z)  # normed integral, no front factors

        front_factor = 3./2.*SOL/(self.fiber.rk/self.fiber.rho)**2
        return coupling_matrix*front_factor