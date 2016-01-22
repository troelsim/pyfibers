from pyfibers.plotter import pgfplotter
from pyfibers.shortcuts import lkfiber as fiber
from pyfibers.coupling.base import get_coupling_matrix
from pyfibers.modes.lekien import LeKienGuidedMode
import numpy as np

N =4  # number of atoms
with pgfplotter('time_evolution') as plt:
    r = 1.1 # Distance from fiber
    g0 = 5.47067930714 # One-atom decay rate

    points = 201

    mode = LeKienGuidedMode(fiber, pol=1, f=1)
    rb = mode.rb

    dz = np.pi/mode.rb*1.1
    zs = [i*dz for i in range(N)]
    phis = [0]*N

    coupling_matrix = get_coupling_matrix(fiber, r, zs, phis)
    evs, evecs = np.linalg.eigh(coupling_matrix)
    expansion = np.linalg.inv(evecs)[:, 6]

    times = np.linspace(0,3, points)/g0
    pops = np.zeros((3*N, points), dtype='complex')
    for i, time in enumerate(times):
        pops[:, i] = np.dot(evecs, expansion*np.exp(-evs*time))

    for state in range(3*N):
        plt.plot(times*g0, abs(pops[state, :])**2, '-k')
    plt.plot(times*g0, sum(abs(pops[s, :])**2 for s in range(3*N)), ':k')
    plt.title("Time evolution of a single-atom excitation")
    plt.ylabel("Population $|e_j^\\nu(t)|^2$")
    plt.xlabel("Time $t/\Gamma_0^{-1}$")
    plt.tight_layout(pad=0)

