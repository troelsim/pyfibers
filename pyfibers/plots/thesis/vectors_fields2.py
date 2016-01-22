from pyfibers.plotter import pgfplotter
from pyfibers.modes.lekien import LeKienGuidedMode
from pyfibers.fibers import LeKienFiber
import numpy as np


@np.vectorize
def rect(r, phi):
    return r*np.cos(phi), r*np.sin(phi)


@np.vectorize
def rect_vec(phi, er, ephi):
    return er*np.cos(phi) - ephi*np.sin(phi), er*np.sin(phi) + ephi*np.cos(phi)


@np.vectorize
def pol(x, y):
    return np.sqrt(x**2 + y**2), np.arctan2(y, x)


with pgfplotter('vector_field_guided') as plt:
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

    plt.figure(figsize=(483.0/72, 624.5/72))
    plt.title('LALA')
    for j in [1,2]:
        mode = LeKienGuidedMode(fiber, 1 if j==1 else -1, 1)
        for i in range(3):
            phase = np.exp(-.25j*np.pi*i)

            e_r_grid = np.real(phase*mode.e_r(Rs, Phis, 0))
            e_phi_grid = np.real(phase*mode.e_phi(Rs, Phis, 0))

            e_x_grid, e_y_grid = rect_vec(Phis, e_r_grid, e_phi_grid)
            e_x_grid = e_x_grid/np.sqrt(mode.closed_norm)
            e_y_grid = e_y_grid/np.sqrt(mode.closed_norm)

            plt.subplot(3, 2, 2*i + j)
            plt.axis('equal')
            plt.quiver(X, Y, e_x_grid, e_y_grid, scale=5)

            if j == 1:
                plt.ylabel({
                    0:  '$\\phi=0$',
                    1:  '$\\phi = \\frac{\\pi}{4}',
                    2:  '$\\phi = \\frac{\\pi}{2}'
                }[i], rotation=0)

            if i == 0:
                plt.title('$p = %s$' % ('+1' if j == 1 else '-1'))

            ax = plt.gca()
            circle = plt.Circle((0, 0), 1, color=None, edgecolor='k', fill=False)
            ax.add_artist(circle)

    plt.tight_layout(pad=0)

with pgfplotter('vector_field_guided_axial') as plt:
    fiber = LeKienFiber(1.45, 1, 1, 1.55)
    mode = LeKienGuidedMode(fiber, -1, 1)

    N = 13
    r = 3
    z = 6

    phase = 1.0+0j

    r_range = np.linspace(0, r, N)
    z_range = np.linspace(0, z, N)

    X, Z = np.meshgrid(r_range, z_range)

    for j in [1, 2]:
        mode = LeKienGuidedMode(fiber, 1, 1 if j == 1 else -1)
        print mode.rb/np.pi
        for i in range(3):
            ax = plt.subplot(2, 3, i+1+3*(j-1))

            phase = np.exp(-1j*np.pi/4*i)

            e_x_grid = np.real(phase*mode.e_r(X, 0, Z))
            e_z_grid = np.real(phase*mode.e_z(X, 0, Z))

            e_x_grid = e_x_grid/np.sqrt(mode.closed_norm)
            e_z_grid = e_z_grid/np.sqrt(mode.closed_norm)

            intensity_grid = np.abs(e_x_grid)**2 + np.abs(e_z_grid)**2
            intensity_grid = intensity_grid/np.max(intensity_grid)

            plt.quiver(X, Z, e_x_grid, e_z_grid, scale=3.5, width=0.006)
            plt.plot([1, 1], [0, z], 'k')
#            plt.axis('equal')

            plt.xlim([0, r])
            plt.ylim([0, z])
            plt.gca().set_aspect('equal', adjustable='box')

            if i == 0:
                plt.ylabel("f=%s" % ('+1' if j==1 else '-1'), rotation=0)
                ax.yaxis.set_label_coords(-0.25, 0.5)

            if j == 1:
                plt.title({
                    0:  '$\\phi=0$',
                    1:  '$\\phi = \\frac{\\pi}{4}',
                    2:  '$\\phi = \\frac{\\pi}{2}'
                }[i], y=1.035)

    plt.tight_layout(pad=0)
