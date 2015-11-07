from pyfibers.modes.lekien import LeKienGuidedMode
from pyfibers.fibers import LeKienFiber
import numpy as np
from matplotlib import pyplot as plt

V = 2.0
for V in np.linspace(0,4,101):
    fiber = LeKienFiber(1.45, 1, 1, V)

    Ws = np.linspace(0, V, 201)
    vals = [LeKienGuidedMode.eigenvalue_equation(fiber, W) for W in Ws]

    plt.plot(Ws, vals)
    plt.interactive(False)
    plt.title('V=%f' % V)
    plt.ylim([-1, 1])
    plt.show()