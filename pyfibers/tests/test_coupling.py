from unittest import TestCase
from pyfibers.modes.lekien import LeKienRadMode
from pyfibers.fibers import LeKienFiber
from pyfibers.atoms import SimpleAtom
from pyfibers.coupling.base import RadCouplingTensor
import numpy as np


class TestCoupling(TestCase):
    def setUp(self):
        self.fiber = LeKienFiber(1.45, 1, 1, 1)
        self.atom_j = SimpleAtom(1.2, 0, 1.5)
        self.atom_k = SimpleAtom(1.2, 1, 0)

    def test_rad_coupling(self):
        coupling = RadCouplingTensor(self.fiber)

        print coupling.get_coupling(self.atom_j, self.atom_k, 0, 1)
        print coupling.calls