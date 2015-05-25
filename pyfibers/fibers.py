from scipy.optimize import brentq
from .modes.lekien import LeKienRadMode, LeKienGuidedMode
import numpy as np
from .constants import *


class Fiber(object):
    rad_mode_class = None
    guided_mode_class = None

    def __init__(self, n, nc, rho, V):
        self.n = n
        self.nc = nc
        self.rho = rho
        self.V = V

    def get_V(self, omega):
        return self.rho*omega*np.sqrt(self.n**2-self.nc**2)/EPSMMU

    def discrete_rad_modes(self, U):
        for mode in self.rad_mode_class.discrete_modes(n=self.n, nc=self.nc, rho=self.rho, V=self.U):
            yield mode

    def discrete_guided_modes(self):
        for mode in self.guided_mode_class.discrete_modes(self):
            yield mode


class LeKienFiber(Fiber):
    rad_mode_class = LeKienRadMode
    guided_mode_class = LeKienGuidedMode
    min_W = 1e-8

    def getUW(self):
        upper_limit  = self.V-self.min_W
        lower_limit = self.min_W

        result, r = brentq(lambda W: self.guided_mode_class.eigenvalue_equation(self, W), lower_limit, upper_limit, full_output=True)

        if r.converged:
            return np.sqrt(self.V**2-result**2), result
        else:
            return None, None
