from pyfibers.fibers import LeKienFiber

def mkradmode(nu, U, pol=1, f=1):
    from pyfibers.fibers import LeKienFiber
    from pyfibers.modes.lekien import LeKienRadMode
    fiber = LeKienFiber(1.45, 1, 1, 1.55)
    return LeKienRadMode(fiber, nu, U, pol, f)


lkfiber = LeKienFiber(1.45, 1, 1, 1.55)
