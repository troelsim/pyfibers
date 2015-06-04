from matplotlib import pyplot as plt, animation
from pyfibers.modes.lekien import LeKienRadMode
from pyfibers.fibers import LeKienFiber
import numpy as np
from copy import deepcopy

radii = np.linspace(0.01, 10, 201)

fiber = LeKienFiber(1.45, 1, 1, 0.77)
#mode = LeKienRadMode(-1, 1.45, 1., 0.77, 1.02)
mode = LeKienRadMode(fiber, 1, 1.02)

e_phis = mode.e_phi(radii, 0, 0)
e_rs = mode.e_r(radii, 0, 0)
e_zs = mode.e_z(radii, 0, 0)

zer = np.linspace(0.01,10,101)
Rer = np.linspace(0.01,2,301)

e_rer = mode.e_r(Rer, 0, 0)
e_phier = mode.e_phi(Rer, 0, 0)
e_zer = mode.e_z(Rer, 0, 0)

plt.figure(1)
plt.plot(Rer, np.real(e_rer))
plt.plot(Rer, np.imag(e_rer))
plt.plot(Rer, 1.45**2*np.real(e_rer))
plt.plot(Rer, 1.45**2*np.imag(e_rer))
plt.figure(2)
plt.plot(Rer, np.real(e_phier))
plt.plot(Rer, np.imag(e_phier))
plt.figure(3)
plt.plot(Rer, np.real(e_zer))
plt.plot(Rer, np.imag(e_zer))

fig, ax = plt.subplots(1, 1)



# Plot of "top view" of fiber vector field
Rs, zs = np.meshgrid(*((np.linspace(0.01,20,21),)*2))

print "rb = ", mode.rb

mode2 = deepcopy(mode)
mode2.nu=mode.nu+1

e_rs = mode.e_r(Rs, 0, zs)+1j*mode2.e_r(Rs, 0, zs)
e_zs = mode.e_z(Rs, 0, zs)+1j*mode2.e_z(Rs, 0, zs)

strength = np.sqrt(abs(e_rs)**2+abs(e_zs)**2)
field = ax.quiver(Rs, zs, e_rs, e_zs, strength, cmap=plt.get_cmap('Blues'))

delta=0.1


def animate(i, field, *args):
    U = np.real(e_rs*np.exp(1j*delta*i))
    V = np.real(e_zs*np.exp(1j*delta*i))
    C = np.sqrt(abs(U)**2+abs(V)**2)

    field.set_UVC(e_rs*np.exp(-1j*delta*i), e_zs*np.exp(-1j*delta*i), C)

anim = animation.FuncAnimation(fig, animate, fargs=(field, Rs, zs), interval=10, blit=False)

plt.show()