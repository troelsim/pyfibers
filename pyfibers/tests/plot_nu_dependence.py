from matplotlib import pyplot as plt, animation
from pyfibers.modes.lekien import LeKienRadMode
from pyfibers.atoms import SimpleAtom
from pyfibers.fibers import LeKienFiber
import numpy as np

radii = np.linspace(0,10,201)

# CW polarization
for nu in range(-2,3):
    mode = LeKienRadMode(nu, 1.45, 1, 2.0+0j, 2.2+0j, pol=1)
    plt.plot(radii, np.real(mode.e_r(radii, 0, 0)), '', label="nu=%d, real part" % nu)
    plt.plot(radii, np.imag(mode.e_r(radii, 0, 0)), '', label="nu=%d, real part" % nu)

plt.legend()
plt.show()

# Animate radial dependence
fiber = LeKienFiber(1, 1.45, 1, 1)
mode = LeKienRadMode(fiber, 1, 1.1, pol=1, f=1)

e_r = mode.e_r(radii, 0, 0)

fig, ax = plt.subplots(1, 1)

graph, = ax.plot(radii, np.real(e_r))

delta = 0.1

#def animate(i):
#    graph.set_ydata(np.real(e_r*np.exp(1j*delta*i)))

#anim = animation.FuncAnimation(fig, animate, interval=10, blit=False)

plt.show()
