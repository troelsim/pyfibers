from pyfibers.fibers import LeKienFiber
from pyfibers.modes.lekien import LeKienRadMode
from ffibers.fibers import fibers
from pyfibers.atoms import SimpleAtom
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import norm
rho = 200.e-9
wavelength = 854e-9
k = 2*np.pi/wavelength

#k = k*rho
#rho = 1.0

n = 1.45
nc = 1.0
V0 = rho*k*np.sqrt(n**2-nc**2)

fiber = LeKienFiber(n, nc, rho, rho*k*np.sqrt(n**2-nc**2))

#fiber = LeKienFiber(1.45, 1, 1, 1)
fibers.init_fiber(1.45, 1, 1, 1)

ms = range(-5,6)
Rs = np.linspace(1,10,201)
es = []
norms = []
es_f = []
norms_f = []
sum_es = Rs*0.0
for m in ms:
    for p in [-1,1]:
        for f in [-1,1]:
            mode = LeKienRadMode(fiber, m, V0*1.01, f=f, pol=p)
            fibers.init_mode(m, 1.3, f, p)

            plt.figure(4)
#            sum_es += abs(mode.e_r(Rs,0,0))**2/mode.norm()
            sum_es += np.array([norm(mode.e_vector(SimpleAtom(R,0,0)))**2 for R in Rs])/mode.norm()
            plt.plot(Rs, np.imag(mode.e_r(Rs,0,0))/np.sqrt(mode.norm()), '', label="%d" % m)
            #plt.plot(Rs, np.vectorize(fibers.e_z)(Rs,0,0))
    es.append(abs(mode.e_z(1,0,0)))
    norms.append(np.sqrt(mode.norm()))

    es_f.append(abs(fibers.e_z(1,0,0)))
    norms_f.append(np.sqrt(abs(fibers.norm())))

plt.legend()
plt.figure(5)
plt.plot(Rs, sum_es)
plt.plot([1],[0])
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

