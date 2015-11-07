from pyfibers.plotter import pgfplotter
from ffibers.fibers import fibers
import numpy as np

fibers.init_fiber(1.45, 1, 1, 1.55)
fibers.init_mode(1, 2.1, 1, 1)

# Plot rad modes
with pgfplotter('rad_e_pol') as plt:
    Rs = np.linspace(0, 8, 301)
    ezsp = np.vectorize(fibers.e_z)(Rs, 0, 0)
    ersp = np.imag(np.vectorize(fibers.e_r)(Rs, 0, 0))
    ephisp = np.vectorize(fibers.e_phi)(Rs, 0, 0)

    fibers.init_mode(-1, 2.1, -1, 1)
    ezsm = np.vectorize(fibers.e_z)(Rs, 0, 0)
    ersm = np.imag(np.vectorize(fibers.e_r)(Rs, 0, 0))
    ephism = np.vectorize(fibers.e_phi)(Rs, 0, 0)

    ylim = (-0.6, 0.6)

    plt.subplot(3,1,1)

    plt.title("Radiative mode, $V=1.55$, $U=2.10$, $m=p=\pm 1$")
    plt.plot(Rs, ersp, 'k', Rs, ersm, ':k')
    plt.ylim(ylim)
    plt.ylabel('$-i e_r$')

    plt.subplot(3,1,2)
    plt.plot(Rs, ephisp, 'k', Rs, ephism, ':k')
    plt.ylim(ylim)
    plt.ylabel('$e_\phi$')

    plt.subplot(3,1,3)
    plt.plot(Rs, ezsp, 'k', Rs, ezsp, ':k')
    plt.ylim(ylim)

    plt.xlabel("$R$")
    plt.ylabel("$e_z$")

with pgfplotter('guid_e_pol') as plt:
    from pyfibers.modes.lekien import LeKienGuidedMode
    from pyfibers.fibers import LeKienFiber
    fiber = LeKienFiber(1.45, 1, 1, 1.55)
    mp = LeKienGuidedMode(fiber, 1, 1)
    mm = LeKienGuidedMode(fiber, -1, 1)
    Rs = np.linspace(0, 8, 301)
    ylim = (-2.5, 2.5)

    plt.subplot(3,1,1)
    plt.title("Guided mode, $V=1.55$, $p=\pm 1$")


    plt.plot(Rs, np.imag(mp.e_r(Rs,0,0)), 'k')
    plt.ylabel('$-i e_r$')
    plt.ylim(ylim)

    plt.subplot(3,1,2)
    plt.plot(Rs, mp.e_phi(Rs, 0, 0), 'k')
    plt.plot(Rs, mm.e_phi(Rs, 0, 0), '--k')
    plt.ylabel('$e_\phi$')
    plt.ylim(ylim)

    plt.subplot(3,1,3)
    plt.plot(Rs, mp.e_z(Rs, 0, 0), 'k')
    plt.ylabel('$e_z$')
    plt.ylim(ylim)

    plt.xlabel('$R=r/\\rho$')


Vs = [1.2, 1.55, 2.0, 4.0]

with pgfplotter('power_contour') as plt:
    from pyfibers.modes.lekien import LeKienGuidedMode
    from pyfibers.fibers import LeKienFiber
    points = 301
    x = np.linspace(-4, 4, points)
    y = np.linspace(-4, 4, points)[:, np.newaxis]
    X, Y = np.meshgrid(x, y)
    plt.figure(figsize=(383.0/72, 383./72))
    R = np.sqrt(X**2+Y**2)
    factor = np.vectorize(lambda r: 1.45**2 if r<1 else 1)(R)
    for i, V in enumerate(Vs):
        plt.subplot(2, 2, i+1)
        fiber = LeKienFiber(1.45, 1, 1, V)
        mp = LeKienGuidedMode(fiber, 1, 1)
        e_rs = mp.e_r(R, 0, 0)*factor
        e_phis = mp.e_phi(R, 0, 0)*factor
        e_zs = mp.e_z(R, 0, 0)*factor
        e = np.sqrt(abs(e_rs)**2 + abs(e_phis)**2 + abs(e_zs)**2)/np.sqrt(mp.closed_norm)
        print np.max(e)
        plt.contourf(X, Y, e**2/(1  0*np.average(e**2)), 100, aspect='equal', cmap=plt.cm.gray_r, vmin=0, vmax=1.45)
        plt.annotate('$V = %.2f$' % V, xy=(-3, 3))
        plt.axis('equal')
