from ffibers.fibers import fibers
import numpy as np
from matplotlib import pyplot as plt

fibers.init_fiber(1.45, 1, 1, 1)
fibers.init_mode(1, 1.3, 1, 1)

Rs = np.linspace(0.1,5,6)

betas = np.linspace(0., 1./np.sqrt(1.45**2-1)-0.01,101)

plt.plot(betas, np.vectorize(fibers.pintegrand)(betas, 0, 0))
plt.show()