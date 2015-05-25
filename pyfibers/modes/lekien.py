from .base import FiberMode
from .util import hankel, hankelp
from scipy.special import jv, jvp, kv, kvp
from scipy.optimize import brentq
from numpy import sqrt, pi, exp

from pyfibers.constants import *


class LeKienRadMode(FiberMode):
    _rho = 1

    @classmethod
    def discrete_modes(cls, m_max=10, n=1.45, nc=1.0, V=2, U=3):
        """
        Iterate over all descrete modal combinations, yielding modes
        :param m_max: max value of the azimuthal wavenumber
        """
        for nu in range(-m_max,m_max+1):
            for pol in [+1, -1]:
                yield cls(nu, n, nc, V, U, pol)


    def __init__(self, nu, n, nc, V, U, pol=1):
        self.nu = nu
        self.pol = pol
        self.n = n
        self.nc = nc
        self.V = V
        Umax = self.rk*self.n
        Umin = self.V
        if U < Umin or U > Umax:
            raise ValueError("U is %f, but must be between %f and %f" % (U, Umin, Umax))
        self.U = U

    def e_r_core(self, R, phi, z):
        return 1j/self.U**2*(
            self.rb*self.U*self.A*jvp(self.nu, self.U*R, 1)
            + 1j*self.nu*self.rk/(R*EPSMMU)*self.B*jv(self.nu, self.U*R)
        )

    def e_r_cladding(self, R, phi, z):
        return 1j/self.Q**2 * \
            sum([
                self.rb*self.Q*self.Cj(j)*hankelp(j, self.nu, self.Q*R) + 1j*self.nu*self.rk/(R*EPSMMU)*self.Dj(j)*hankel(j, self.nu, self.Q*R)
                for j in [1, 2]
            ])

    def e_phi_core(self, R, phi, z):
        return 1j/self.U**2*(
            1j*self.nu*self.rb/R*self.A*jv(self.nu, self.U*R) - self.U*self.rk/(EPSMMU)*self.B*jvp(self.nu, self.U*R)
        )

    def e_phi_cladding(self, R, phi, z):
        return 1j/self.Q**2*sum([
            1j*self.nu*self.rb/R*self.Cj(j)*hankel(j, self.nu, self.Q*R)
            - self.Q*self.rk/EPSMMU*self.Dj(j)*hankelp(j, self.nu, self.Q*R)
            for j in [1,2]
        ])

    def e_z_core(self, R, phi, z):
        return self.A*jv(self.nu, self.U*R)

    def e_z_cladding(self, R, phi, z):
        return sum([
            self.Cj(j)*hankel(j, self.nu, self.Q*R)
            for j in [1, 2]
        ])

    def implicit(self, R, phi, z):
        return exp(1j*self.rb/self.rho*z+1j*self.nu*phi)

    def norm(self):
        return 2*pi*self.rk*self.rho/self.Q**2/EPSMMU*(self.nc**2*abs(self.Cj(1))**2+MU0/EPS0*abs(self.Dj(1))**2)


    @property
    def A(self):
         return 1

    @property
    def B(self):
        return 1j*self.pol*self.eta(2)*self.A

    def eta(self, j):
        return EPS0*SOL*sqrt((self.nc*abs(self.Vj(j))**2 + abs(self.Lj(j))**2)/(abs(self.Vj(j))**2 + self.nc*abs(self.Mj(j))**2))

    @property
    def Q(self):
        return sqrt(self.U**2-self.V**2)

    @property
    def rb(self):
        Q = self.Q
        U = self.U
        return sqrt((U**2*self.nc**2 - Q**2*self.n**2)/(self.n**2-self.nc**2))

    @property
    def rk(self):
        return self.V/sqrt(self.n**2-self.nc**2)

    @property
    def rho(self):
        return self._rho or 1  # TODO what shall we do about this?

    def Vj(self, j):
        return self.nu*self.rho*self.rb*self.rk / \
            (self.Q**2*self.U**2)*(self.nc**2-self.n**2)*jv(self.nu, self.U)*hankel(j, self.nu, self.Q).conjugate()

    def Mj(self, j):
        return self.rho/self.U*jvp(self.nu, self.U, 1)*hankel(j, self.nu, self.Q).conjugate() - \
            self.rho/self.Q*jv(self.nu, self.U)*hankelp(j, self.nu, self.Q).conjugate()

    def Lj(self, j):
        return self.rho*self.n**2/self.U*jvp(self.nu, self.U, 1)*hankel(j, self.nu, self.Q).conjugate() - \
            self.rho/self.Q*self.nc**2*jv(self.nu, self.U)*hankelp(j, self.nu, self.Q).conjugate()

    def Cj(self, j):
        return (-1)**j * 1j*pi*self.Q**2/(4*self.rho*self.nc**2)*(self.A*self.Lj(j) + 1j*MU0*SOL*self.B*self.Vj(j))

    def Dj(self, j):
        return (-1)**(j-1) * 1j*pi*self.Q**2/(4*self.rho)*(1j*EPS0*SOL*self.A*self.Vj(j) - self.B*self.Mj(j))


class LeKienGuidedMode(FiberMode):

    min_W = 1e-10

    def __init__(self, fiber, pol, f):
        self.fiber = fiber
        self.n = fiber.n
        self.nc = fiber.nc
        self.V = fiber.V
        self.U, self.W = self.getUW(fiber)
        self.pol = pol
        self.f = f

    @classmethod
    def discrete_modes(cls, fiber, **kwargs):
        for pol in [1,-1]:
            for f in [1,-1]:
                yield cls(fiber, pol, f)

    @classmethod
    def eigenvalue_equation(cls, fiber, W):
        U = sqrt(fiber.V**2 - W**2)
        j0 = jv(0, U)
        j1 = jv(1, U)
        k0 = kv(0, W)
        k1 = kv(1, W)
        k1p = kvp(1, W, 1)
        return (
            j0/(U*j1) + (fiber.n**2 + fiber.nc**2)/(2*fiber.n**2) * k1p/(W*k1) - 1/U**2 +
            sqrt(
                ((fiber.n**2 - fiber.nc**2)/(2*fiber.n**2) * k1p/(W*k1))**2
                + (cls._rboverrk(fiber, U, fiber.V)**2/fiber.n**2) * (1/W**2 + 1/U**2)**2
            )
        )

    @staticmethod
    def _rboverrk(fiber, U, V):
        return sqrt(
            (fiber.n**2*V**2 - (fiber.n**2-fiber.nc**2)*U**2) /
            V**2
        )

    @property
    def rb(self):
        return sqrt(
            (self.fiber.n**2 * self.W**2 + self.fiber.nc**2 * self.U**2) /
            (self.fiber.n**2 - self.fiber.nc**2)
        )

    @classmethod
    def getUW(cls, fiber):
        upper_limit  = fiber.V-cls.min_W
        lower_limit = cls.min_W

        result, r = brentq(lambda W: cls.eigenvalue_equation(fiber, W), lower_limit, upper_limit, full_output=True)

        if r.converged:
            return sqrt(fiber.V**2-result**2), result
        else:
            return None, None

    @property
    def s(self):
        U = self.U
        W = self.W
        return (1/U**2 + 1/W**2)/(jvp(1, U, 1)/(U*jv(1,U)) + kvp(1, W, 1)/(W*kv(1, W)))

    def e_r_core(self, R, phi, z):
        U = self.U
        W = self.W
        return 1j*self.W/self.U*kv(1, W)/jv(1, U)*(
            (1-self.s)*jv(0, U*R) - (1+self.s)*jv(2, U*R)
        )

    def e_r_cladding(self, R, phi, z):
        U = self.U
        W = self.W
        return 1j*(
            (1-self.s)*kv(0, W*R) + (1+self.s)*kv(2, W*R)
        )

    def e_phi_core(self, R, phi, z):
        U = self.U
        W = self.W
        return - W/U * kv(1, W)/jv(1, U) * (
            (1-self.s)*jv(0, U*R) + (1+self.s)*jv(2, U*R)
        )

    def e_phi_cladding(self, R, phi, z):
        U = self.U
        W = self.W
        return - (
            (1-self.s)*kv(0, W*R) - (1+self.s)*kv(2, W*R)
        )

    def e_z_core(self, R, phi, z):
        return 2*self.W/self.rb * kv(1, self.W)/jv(1, self.U)*jv(1, self.U*R)

    def e_z_cladding(self, R, phi, z):
        return 2*self.W/self.rb * kv(1, self.W*R)

    def implicit(self, R, phi, z):
        return 1+0j