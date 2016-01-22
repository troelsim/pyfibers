from pyfibers.plotter import pgfplotter
from pyfibers.shortcuts import lkfiber as fiber
from pyfibers.coupling.base import get_coupling_matrix
from pyfibers.modes.lekien import LeKienGuidedMode
import numpy as np
#from matplotlib import pyplot as plt

N =4  # number of atoms
with pgfplotter('guided_rates') as plt:
    plt.figure(figsize=(383.0/72, 383./72))
    for pic, N in enumerate([2, 4, 10, 20]):
        r = 1.1 # Distance from fiber
        g0 = 5.47067930714 # One-atom decay rate

        points = 201

        mode = LeKienGuidedMode(fiber, pol=1, f=1)
        rb = mode.rb

        dzs = np.linspace(0, 2, points)*np.pi/rb
        evs_total = np.zeros((3*N, points))  # To hold the results

        for j, dz in enumerate(dzs):  # For each atomic separation calculate the eigenvalues
            zs = [i*dz for i in range(N)]
            phis = [0]*N

            coupling_matrix = get_coupling_matrix(fiber, r, zs, phis)
            evs, evecs = np.linalg.eigh(coupling_matrix)

            evs_total[:, j] = evs/g0

        plt.subplot(2,2,pic+1)
        for i in range(3*N):
            plt.plot(dzs*rb/np.pi, evs_total[i, :], 'k-')
        plt.plot(dzs*rb/(np.pi), np.sum(evs_total, axis=0), 'k:', label="SUM")
        m = np.sum(evs_total[:,0])
        plt.xlabel("$\Delta z  \\beta / \pi$")
        plt.ylabel("$\Gamma/\Gamma_0$")
        plt.title("$N=%d$" % N)
        plt.ylim([0, m*1.1])
    plt.tight_layout(pad=0)
