from pyfibers.plotter import pgfplotter
from pyfibers.modes.lekien import LeKienGuidedMode
from pyfibers.fibers import LeKienFiber
import numpy as np
from scipy.stats import gamma
from scipy.optimize import curve_fit


with pgfplotter('eta_guided') as plt:
    Vs = np.linspace(0.1, 10, 101)
    etas = []
    for V in Vs:
        print V
        fiber = LeKienFiber(1.45, 1, 1, V)
        try:
            mode = LeKienGuidedMode(fiber, 1, 1)
            print "Found U=%f, W=%f" % (mode.U, mode.W)
        except LeKienGuidedMode.NoSolutionsFoundException:
            etas.append(0)
            print "no luck"
            continue
        etas.append(mode.eta)

    plt.plot(Vs, etas, 'k')
    plt.title("Fraction of power transported inside fiber core")
    plt.xlabel("$V$")
    plt.ylabel("$P_1/(P_1+P_2)$")
