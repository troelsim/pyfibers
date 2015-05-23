from .base import FiberMode
from .util import hankel, hankelp
from scipy.special import jv, jvp
from numpy import sqrt, pi, exp


EPS0 = 8.854187817e-12
MU0 = 4*pi*1e-7
SOL = 299792458  # Speed of light
EPSMMU = sqrt(EPS0/MU0)

class LeKienRadMode(FiberMode):
    _rho = 1

    def __init__(self, nu, n, nc, V, U, pol=1):
        self.nu = nu
        self.pol = pol
        self.n = n
        self.nc = nc
        self.V = V
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
        Q = self.Q
        U = self.U
        return sqrt((U**2-Q**2)/(self.n**2-self.nc**2))

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




