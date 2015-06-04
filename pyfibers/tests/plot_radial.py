import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import LeKienRadMode
from pyfibers.modes.util import hankel

fiber = LeKienFiber(1.45, 1, 1, 1)
#m = LeKienRadMode(-5, 1.45, 1, 0.77, 1.06, pol=1)
m = LeKienRadMode(fiber, 1, 1.06, pol=1)

Rs = np.linspace(0.01, 20, 201)
#e_zs = hankel(2, 8, Rs)
e_zs = m.e_r(Rs, 0, 0)
m.pol = -1
e_zs += 1j*m.e_r(Rs, 0, 0)


fig = plt.figure()
ax = plt.axes()
line1, = ax.plot(Rs, np.real(e_zs))
line2, = ax.plot(Rs, np.imag(e_zs))
line3, = ax.plot(Rs, np.abs(e_zs))

#plt.show()


def animate(i):
    y = e_zs*np.exp(1j*0.09*i)
    line1.set_data(Rs, np.real(y))
    line2.set_data(Rs, np.imag(y))
    line3.set_data(Rs, np.abs(y)**2)
    return line1, line2

anim = animation.FuncAnimation(fig, animate, interval=1)

plt.show()
