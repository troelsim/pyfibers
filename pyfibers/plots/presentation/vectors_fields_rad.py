#from pyfibers.plotter import pgfplotter
from pyfibers.modes.lekien import LeKienGuidedMode
from pyfibers.fibers import LeKienFiber
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation


@np.vectorize
def rect(r, phi):
    return r*np.cos(phi), r*np.sin(phi)


@np.vectorize
def rect_vec(phi, er, ephi):
    return er*np.cos(phi) - ephi*np.sin(phi), er*np.sin(phi) + ephi*np.cos(phi)


@np.vectorize
def pol(x, y):
    return np.sqrt(x**2 + y**2), np.arctan2(y, x)


fiber = LeKienFiber(1.45, 1, 1, 1.55)
mode = LeKienGuidedMode(fiber, -1, 1)

N = 15
r = 2

phase = np.exp(0.0j)

x_range = np.linspace(-r, r, N)
y_range = np.linspace(-r, r, N)

X, Y = np.meshgrid(x_range, y_range)
Rs, Phis = pol(X, Y)

e_r_grid = np.real(phase*mode.e_r(Rs, Phis, 0))
e_phi_grid = np.real(phase*mode.e_phi(Rs, Phis, 0))

e_x_grid, e_y_grid = rect_vec(Phis, e_r_grid, e_phi_grid)
e_x_grid = e_x_grid/np.sqrt(mode.closed_norm)
e_y_grid = e_y_grid/np.sqrt(mode.closed_norm)

fig = plt.figure(figsize=(483.0/72, 483.0/72))

fpc = 50  # Frames per cycle

for j in [1,2]:
    mode = LeKienGuidedMode(fiber, 1 if j==1 else -1, 1)
    def update_quiver(frame, quiver):
        phase = np.exp(2.j*np.pi*frame/fpc)
        e_r_grid = np.real(phase*mode.e_r(Rs, Phis, 0))
        e_phi_grid = np.real(phase*mode.e_phi(Rs, Phis, 0))

        e_x_grid, e_y_grid = rect_vec(Phis, e_r_grid, e_phi_grid)
        e_x_grid = e_x_grid/np.sqrt(mode.closed_norm)
        e_y_grid = e_y_grid/np.sqrt(mode.closed_norm)
        quiver.set_UVC(e_x_grid, e_y_grid)

    fig = plt.figure(figsize=(483.0/72, 483.0/72))
    ax = plt.gca()

    circle = plt.Circle((0, 0), 1, color=None, edgecolor='k', fill=False)
    ax.add_artist(circle)

    Q = ax.quiver(X, Y, e_x_grid, e_y_grid, scale=5)
    anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q,), interval=40, blit=False)
    plt.axis('equal')
    anim.save('guided-cross-section-pol%d.mp4' % j, fps=25, bitrate=1024)
    #plt.show()

import sys
#sys.exit()

fiber = LeKienFiber(1.45, 1, 1, 1.55)
mode = LeKienGuidedMode(fiber, -1, 1)

N = 13
r = 3
z = 6

phase = 1.0+0j

r_range = np.linspace(0, r, N)
z_range = np.linspace(0, z, N)

X, Z = np.meshgrid(r_range, z_range)

e_x_grid = np.real(phase*mode.e_r(X, 0, Z))
e_z_grid = np.real(phase*mode.e_z(X, 0, Z))

e_x_grid = e_x_grid/np.sqrt(mode.closed_norm)
e_z_grid = e_z_grid/np.sqrt(mode.closed_norm)

for f in [1, -1]:
    mode = LeKienGuidedMode(fiber, -1, f)

    def update_quiver(frame, quiver):
        phase = np.exp(2.j*np.pi*frame/fpc)
        e_x_grid = np.real(phase*mode.e_r(X, 0, Z))
        e_z_grid = np.real(phase*mode.e_z(X, 0, Z))

        e_x_grid = e_x_grid/np.sqrt(mode.closed_norm)
        e_z_grid = e_z_grid/np.sqrt(mode.closed_norm)
        quiver.set_UVC(e_x_grid, e_z_grid)

    fig = plt.figure(figsize=(483.0/72, 2*483.0/72))
    ax = plt.gca()

    Q = ax.quiver(X, Z, e_x_grid, e_z_grid, scale=2, width=0.006)
    ax.plot([1, 1], [0, z], 'k')
    ax.set_aspect('equal', adjustable='box')
    anim = animation.FuncAnimation(fig, update_quiver, fargs=(Q,), interval=40, blit=False)
    anim.save('guided-axial-f%d.mp4' % f, fps=25, bitrate=1024)
