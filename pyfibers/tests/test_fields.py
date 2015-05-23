from unittest import TestCase
from modes.lekien import LeKienRadMode
import numpy as np

class TestLeKienFields(TestCase):
    def setUp(self):
        V = 2.0+0j
        U = 2.1+0j
        self.nc = 1.0
        self.n = 1.45
        self.mode = LeKienRadMode(1, self.n, self.nc, V, U)

    def test_er(self):
        self.assertAlmostEqual(self.mode.e_r_cladding(1,0,0)*self.nc**2, self.mode.e_r_core(1,0,0)*self.n**2)

    def test_ephi(self):
        self.assertAlmostEqual(self.mode.e_phi_cladding(1,0,0), self.mode.e_phi_core(1,0,0))

    def test_ez(self):
        self.assertAlmostEqual(self.mode.e_z_cladding(1,0,0), self.mode.e_z_core(1,0,0))
