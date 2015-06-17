from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import LeKienRadMode
from ffibers.fibers import fibers
from matplotlib import pyplot as plt
import numpy as np

fiber = LeKienFiber(1.45, 1, 1, 1)
fibers.init_fiber(1.45, 1, 1, 1)

ms = range(-14,15)
Rs = np.linspace(0,40,201)
es = []
norms = []
es_f = []
norms_f = []
sum_es = Rs*0.0
for m in ms:
    for p in [-1,1]:
        for f in [-1,1]:
            mode = LeKienRadMode(fiber, m, 1.3, f=f, pol=p)
            fibers.init_mode(m, 1.3, f, p)

            plt.figure(4)
            sum_es += abs(mode.e_z(Rs,0,0))**2/mode.norm()
            plt.plot(Rs, mode.e_z(Rs,0,0)/np.sqrt(mode.norm()), '', label="%d" % m)
            #plt.plot(Rs, np.vectorize(fibers.e_z)(Rs,0,0))
    es.append(abs(mode.e_z(1,0,0)))
    norms.append(np.sqrt(mode.norm()))

    es_f.append(abs(fibers.e_z(1,0,0)))
    norms_f.append(abs(fibers.norm()))

plt.legend()
plt.figure(5)
plt.plot(Rs, sum_es)
plt.show()


plt.figure(1)
plt.legend()
plt.semilogy()
plt.plot(ms,es)
plt.plot(ms, es_f)
plt.figure(2)
plt.legend()
plt.semilogy()
plt.plot(ms,norms)
plt.plot(ms, norms_f)
plt.figure(3)
plt.legend()
plt.plot(ms,np.array(es)/np.array(norms))
plt.plot(ms,np.array(es_f)/np.array(norms_f))

plt.show()

