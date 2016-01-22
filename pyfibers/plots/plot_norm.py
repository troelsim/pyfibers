from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import LeKienRadMode
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt

n = 101
m= 201
Vs = np.linspace(0.1, 60, n)
Rs = np.linspace(0.1, 50, m)
norms = []
strengths = np.zeros(m)

for V in Vs:
    fiber = LeKienFiber(1.45, 1, 1.0, V)
    mode = LeKienRadMode(fiber, 1, V*1.1, 1, 1)
    norms.append(mode.norm()*mode.Q**2)
    ezs = mode.e_z(Rs, 0, 0)/np.sqrt(mode.norm())
    ers = mode.e_r(Rs, 0, 0)/np.sqrt(mode.norm())
    ephis = mode.e_phi(Rs, 0, 0)/np.sqrt(mode.norm())
    estrengths = np.sqrt(abs(ezs)**2+abs(ers)**2+abs(ephis)**2)/np.sqrt(mode.norm())
    strengths += estrengths
#    plt.plot(Rs, ezs, label='e_z')
#    plt.semilogy()

plt.plot(Rs, strengths)


#plt.plot(Vs, norms)
plt.show()
