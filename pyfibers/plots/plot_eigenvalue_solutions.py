from pyfibers.modes.lekien import LeKienGuidedMode
from pyfibers.fibers import LeKienFiber
import numpy as np
from matplotlib import pyplot as plt


Vs = np.linspace(0, 10, 201)
Uss = np.zeros((20,201))

for j, V in enumerate(Vs):
    fiber = LeKienFiber(1.45, 1, 1, V)
    mesh = np.linspace(0, V, 1001)
    Us=[]
    for i, U_lower in enumerate(mesh):
        try:
            U, W = LeKienGuidedMode.getUW(fiber, lower_limit=U_lower, upper_limit=mesh[i+1])
            if U and not np.isnan(U):
                Us.append(U)
        except ValueError:
            pass
        except IndexError:
            break
        Uss[:len(Us),j] = Us

plt.plot(Vs, Uss[0,:])
plt.plot(Vs, Uss[1,:])
plt.plot(Vs, Uss[2,:])
plt.plot(Vs, Uss[3,:])
plt.show()