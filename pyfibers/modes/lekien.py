from .base import FiberMode
from .util import hankel, hankelp, cache_for_fiber
from scipy.special import jv, jvp, kv, kvp
from scipy.optimize import brentq
from scipy.integrate import quad
from numpy import sqrt, pi, exp, inf, linspace

from pyfibers.constants import *


class LeKienRadMode(FiberMode):

    @classmethod
    def discrete_modes(cls, fiber, m_max=10, U=3):
        """
        Iterate over all descrete modal combinations, yielding modes
        :param m_max: max value of the azimuthal wavenumber
        """
        for nu in range(-m_max,m_max+1):
            for pol in [+1, -1]:
                for f in [+1, -1]:
                    yield cls(fiber, nu, U, pol, f=f)

    def __init__(self, fiber, nu, U, pol=1, f=1):
        self.fiber = fiber
        self.nu = nu
        self.pol = pol
        self.n = fiber.n
        self.nc = fiber.nc
        self.V = fiber.V
        self.rho = fiber.rho
        self.f = f
        Umax = self.rk*self.n
        Umin = self.V
        if U < Umin or U > Umax:
            raise ValueError("U is %f, but must be between %f and %f" % (U, Umin, Umax))
        self.U = U

    @staticmethod
    def get_U(fiber, rb):
        return sqrt(fiber.n**2*fiber.rk**2 - rb**2)


    @staticmethod
    def U_max(fiber):
        return fiber.rk*fiber.n

    @staticmethod
    def U_min(fiber):
        return fiber.V

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
        return exp(1j*self.f*self.rb/self.rho*z+1j*self.nu*phi)

    def norm(self):
        return 8*pi*self.rk*self.rho*SOL/self.Q**2*(self.nc**2*abs(self.Cj(1))**2+MU0/EPS0*abs(self.Dj(1))**2)


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
        return self.f*sqrt((U**2*self.nc**2 - Q**2*self.n**2)/(self.n**2-self.nc**2))

    @property
    def rk(self):
        return self.V/sqrt(self.n**2-self.nc**2)

    @cache_for_fiber
    def Vj(self, j):
        return self.nu*self.rb*self.rk / \
            (self.Q**2*self.U**2)*(self.nc**2-self.n**2)*jv(self.nu, self.U)*hankel(j, self.nu, self.Q).conjugate()

    @cache_for_fiber
    def Mj(self, j):
        return 1./self.U*jvp(self.nu, self.U, 1)*hankel(j, self.nu, self.Q).conjugate() - \
            1./self.Q*jv(self.nu, self.U)*hankelp(j, self.nu, self.Q).conjugate()

    @cache_for_fiber
    def Lj(self, j):
        return 1.*self.n**2/self.U*jvp(self.nu, self.U, 1)*hankel(j, self.nu, self.Q).conjugate() - \
            1./self.Q*self.nc**2*jv(self.nu, self.U)*hankelp(j, self.nu, self.Q).conjugate()

    @cache_for_fiber
    def Cj(self, j):
        return (-1)**j * 1j*pi*self.Q**2/(4*self.nc**2)*(self.A*self.Lj(j) + 1j*MU0*SOL*self.B*self.Vj(j))

    @cache_for_fiber
    def Dj(self, j):
        return (-1)**(j-1) * 1j*pi*self.Q**2/4*(1j*EPS0*SOL*self.A*self.Vj(j) - self.B*self.Mj(j))


class LeKienGuidedMode(FiberMode):

    min_W = 1e-10

    class NoSolutionsFoundException(BaseException):
        def __init__(self, s=""):
            self.message = s

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
    @cache_for_fiber
    def rb(self):
        return sqrt(
            (self.fiber.n**2 * self.W**2 + self.fiber.nc**2 * self.U**2) /
            (self.fiber.n**2 - self.fiber.nc**2)
        )

    @property
    @cache_for_fiber
    def rk(self):
        return self.fiber.V/sqrt(self.fiber.n**2-self.fiber.nc**2)

    @classmethod
    def getUW(cls, fiber):
        """
        We have to solve the eigenvalue equaition, which can have a singularity
        and zero or more solutions.
        As it turns out, the right solution is when the function is increasing.
        The at the singularity s, f goes to infinity for x->s+ and -inf for x->s-, meaning
        an interval (a, b) will contain the solution if f(a)<0 and f(b)>0 and a and b are close
        :param fiber:
        :return:
        """
        L = 2.405
        upper_limit = fiber.V-cls.min_W
        lower_limit = cls.min_W
        if fiber.V > L:
            lower_limit = sqrt(fiber.V**2-L**2)

        test_ws = linspace(lower_limit, upper_limit, 101)

        previous_val = 0
        previous_w = 0

        solutions = []

        for w in test_ws:
            val = cls.eigenvalue_equation(fiber, w)
            if previous_val * val < 0:
                result, r = brentq(lambda W: cls.eigenvalue_equation(fiber, W), previous_w, w, full_output=True)
                unfitness = (abs(val)+abs(previous_val))  # Unlikely to be a root if this is high
                if r.converged:
                    #print "unfitness %f" % unfitness
                    solutions.append((sqrt(fiber.V**2 - result**2), result, unfitness))
                else:
                    raise cls.NoSolutionsFoundException("Brent could not help you")
            previous_val = val
            previous_w = w
        if solutions:
            #print "%d solution(s)" % len(solutions)
            for U, W, unfitness in solutions[::-1]:
                if unfitness < 1:
                    return U, W

        raise cls.NoSolutionsFoundException("No solution found")

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
        return -self.pol * W/U * kv(1, W)/jv(1, U) * (
            (1-self.s)*jv(0, U*R) + (1+self.s)*jv(2, U*R)
        )

    def e_phi_cladding(self, R, phi, z):
        U = self.U
        W = self.W
        return -self.pol * (
            (1-self.s)*kv(0, W*R) - (1+self.s)*kv(2, W*R)
        )

    def e_z_core(self, R, phi, z):
        return 2*self.f*self.W/self.rb * kv(1, self.W)/jv(1, self.U)*jv(1, self.U*R)

    def e_z_cladding(self, R, phi, z):
        return 2*self.f*self.W/self.rb * kv(1, self.W*R)

    def implicit(self, R, phi, z):
        return exp(1j*self.f*self.rb/self.fiber.rho*z+1j*self.pol*phi)

    @cache_for_fiber
    def norm(self):
        result, err = quad(
            lambda R: 2*pi*R*(
                (self.n**2 if R<1 else self.nc**2)
                * sqrt(
                    abs(self.e_z(R,0,0))**2 +
                    abs(self.e_phi(R,0,0))**2 +
                    abs(self.e_z(R,0,0))**2
                )*self.fiber.rho**3
            ),
            0,
            inf
        )
        return result

    @property
    @cache_for_fiber
    def P1(self):
        U = self.U
        W = self.W
        return W**2/U**2 * kv(1, W)**2/jv(1, U)**2 * (
            (1-self.s)**2 * (jv(0, U)**2 + jv(1, U)**2) +
            (1+self.s)**2 * (jv(2, U)**2 - jv(1, U)*jv(3, U)) +
            2*U**2/self.rb**2 * (jv(1, U)**2 - jv(0, U)*jv(2, U))
        )

    @property
    @cache_for_fiber
    def P2(self):
        U = self.U
        W = self.W
        return (1-self.s)**2 * (kv(1, W)**2 - kv(0, W)**2) +\
            (1+self.s)**2 * (kv(1, W)*kv(3, W) - kv(2, W)**2) +\
            2*W**2/self.rb**2*(kv(0, W)*kv(2, W) - kv(1, W)**2)

    @property
    @cache_for_fiber
    def closed_norm(self):
        return 2*pi * self.fiber.rho**2 * (self.n*self.P1 + self.nc*self.P2)

    @property
    @cache_for_fiber
    def betaprime(self):
        return self.fiber.n**2*self.rk/self.rb/SOL * (1-2*self.fiber.delta*(1-self.P1/(self.P1+self.P2)))

    @property
    @cache_for_fiber
    def eta(self):
        return self.P1/(self.P2 + self.P1)
