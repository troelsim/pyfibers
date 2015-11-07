from functools import wraps
from unittest import TestCase
from ffibers.fibers import fibers
from math import sqrt
from scipy.integrate import quad
import numpy as np

def each_mode(func):
    @wraps(func)
    def inner(obj, *args, **kwargs):
        for mode in obj.modes:
            fibers.init_mode(*mode)
            func(obj, *mode)
    return inner

class TestFortranFields(TestCase):
    EPS = 1e-8

    def setUp(self):
        self.n = 1.45
        self.nc = 1.0
        self.phi = 0
        self.z = 0
        fibers.init_fiber(self.n, self.nc, 1, 1)
        U = 1.3
        self.modes = []
        for m in range(-5, 6):
            for pol in [-1, 1]:
                for f in [-1, 1]:
                    self.modes.append((m, U, pol, f))

    @each_mode
    def test_er(self, *args):
        self.assertAlmostEqual(
            fibers.e_r(1-self.EPS, self.phi, self.z) * self.n**2,
            fibers.e_r(1+self.EPS, self.phi, self.z) * self.nc**2
        )

    @each_mode
    def test_ephi(self, *args):
        self.assertAlmostEqual(
            fibers.e_phi(1-self.EPS, self.phi, self.z),
            fibers.e_phi(1+self.EPS, self.phi, self.z)
        )

    @each_mode
    def test_ez(self, *args):
        self.assertAlmostEqual(
            fibers.e_z(1-self.EPS, self.phi, self.z),
            fibers.e_z(1+self.EPS, self.phi, self.z)
        )

    @each_mode
    def test_ml_symmetry_r(self, *args):
        normalization = sqrt(fibers.norm())
        e_r_cladding = fibers.e_r(1+self.EPS, self.phi, self.z)/normalization
        e_r_core = fibers.e_r(1-self.EPS, self.phi, self.z)/normalization
        sign = (-1)**args[0]
        args = (-args[0], args[1], -args[2], args[3])  # Flip m and pol
        fibers.init_mode(*args)
        normalization = sqrt(fibers.norm())
        e_r_cladding_flipped = fibers.e_r(1+self.EPS, self.phi, self.z)/normalization
        e_r_core_flipped = fibers.e_r(1-self.EPS, self.phi, self.z)/normalization
        self.assertAlmostEqual(e_r_cladding, sign*e_r_cladding_flipped)
        self.assertAlmostEqual(e_r_core, sign*e_r_core_flipped)

    @each_mode
    def test_ml_symmetry_phi(self, *args):
        normalization = sqrt(fibers.norm())
        e_phi_cladding = fibers.e_phi(1+self.EPS, self.phi, self.z)/normalization
        e_phi_core = fibers.e_phi(1-self.EPS, self.phi, self.z)/normalization
        sign = (-1)**(args[0]+1)
        args = (-args[0], args[1], -args[2], args[3])  # Flip m and pol
        fibers.init_mode(*args)
        normalization = sqrt(fibers.norm())
        e_phi_cladding_flipped = fibers.e_phi(1+self.EPS, self.phi, self.z)/normalization
        e_phi_core_flipped = fibers.e_phi(1-self.EPS, self.phi, self.z)/normalization
        self.assertAlmostEqual(e_phi_cladding, sign*e_phi_cladding_flipped)
        self.assertAlmostEqual(e_phi_core, sign*e_phi_core_flipped)

    @each_mode
    def test_ml_symmetry_z(self, *args):
        normalization = sqrt(fibers.norm())
        e_z_cladding = fibers.e_z(1+self.EPS, self.phi, self.z)/normalization
        e_z_core = fibers.e_z(1-self.EPS, self.phi, self.z)/normalization
        sign = (-1)**(args[0])
        args = (-args[0], args[1], -args[2], args[3])  # Flip m and pol
        fibers.init_mode(*args)
        normalization = sqrt(fibers.norm())
        e_z_cladding_flipped = fibers.e_z(1+self.EPS, self.phi, self.z)/normalization
        e_z_core_flipped = fibers.e_z(1-self.EPS, self.phi, self.z)/normalization
        self.assertAlmostEqual(e_z_cladding, sign*e_z_cladding_flipped)
        self.assertAlmostEqual(e_z_core, sign*e_z_core_flipped)


    @each_mode
    def test_fl_symmetry_z(self, *args):
        normalization = sqrt(fibers.norm())
        e_z_cladding = fibers.e_z(1+self.EPS, self.phi, self.z)/normalization
        e_z_core = fibers.e_z(1-self.EPS, self.phi, self.z)/normalization
        sign = 1
        args = (args[0], args[1], -args[2], -args[3])  # Flip f and pol
        fibers.init_mode(*args)
        normalization = sqrt(fibers.norm())
        e_z_cladding_flipped = fibers.e_z(1+self.EPS, self.phi, self.z)/normalization
        e_z_core_flipped = fibers.e_z(1-self.EPS, self.phi, self.z)/normalization
        self.assertAlmostEqual(e_z_cladding, sign*e_z_cladding_flipped)
        self.assertAlmostEqual(e_z_core, sign*e_z_core_flipped)

    @each_mode
    def test_fl_symmetry_phi(self, *args):
        normalization = sqrt(fibers.norm())
        e_phi_cladding = fibers.e_phi(1+self.EPS, self.phi, self.z)/normalization
        e_phi_core = fibers.e_phi(1-self.EPS, self.phi, self.z)/normalization
        sign = -1
        args = (args[0], args[1], -args[2], -args[3])  # Flip f and pol
        fibers.init_mode(*args)
        normalization = sqrt(fibers.norm())
        e_phi_cladding_flipped = fibers.e_phi(1+self.EPS, self.phi, self.z)/normalization
        e_phi_core_flipped = fibers.e_phi(1-self.EPS, self.phi, self.z)/normalization
        self.assertAlmostEqual(e_phi_cladding, sign*e_phi_cladding_flipped)
        self.assertAlmostEqual(e_phi_core, sign*e_phi_core_flipped)

    @each_mode
    def test_fl_symmetry_r(self, *args):
        normalization = sqrt(fibers.norm())
        e_r_cladding = fibers.e_r(1+self.EPS, self.phi, self.z)/normalization
        e_r_core = fibers.e_r(1-self.EPS, self.phi, self.z)/normalization
        sign = -1
        args = (args[0], args[1], -args[2], -args[3])  # Flip f and pol
        fibers.init_mode(*args)
        normalization = sqrt(fibers.norm())
        e_r_cladding_flipped = fibers.e_r(1+self.EPS, self.phi, self.z)/normalization
        e_r_core_flipped = fibers.e_r(1-self.EPS, self.phi, self.z)/normalization
        self.assertAlmostEqual(e_r_cladding, sign*e_r_cladding_flipped)
        self.assertAlmostEqual(e_r_core, sign*e_r_core_flipped)

    @each_mode
    def test_norm(self, *args):
        # This makes no sense - the integral will always diverge.
        norm = fibers.norm()

        def integrand(r):
            return 2*np.pi*r*(self.n**2 if r<0 else self.nc**2)*sum([
                abs(getattr(fibers, 'e_'+dim)(r, self.phi, self.z))**2
                for dim in ['r', 'phi', 'z']])

        integral, error = quad(integrand, 0, np.inf, limit=1000)
        print u"I = %e +/- %e" % (integral, error)
        self.assertAlmostEqual(
            norm,
            integral
        )
