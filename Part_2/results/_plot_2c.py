import numpy as np
import matplotlib.pylab as plt

#----------------------------------------
DATA = np.loadtxt("./temp_sweep_128_nearcrit.dat")
nsamples = 100000
n=128

Tcrit = 2 / np.log(1+np.sqrt(2))

T = DATA[:,0]

EN, ENVAR = DATA[:,1], DATA[:,2]
SPIN, SPINVAR = DATA[:,3], DATA[:,4]
SPIN = abs(SPIN)
EN = abs(EN)

#----------------------------------------
EN /= n**2
SPIN/=n**2

ENVAR/=n**4
SPINVAR/=n**4

MS = SPINVAR / T * n**2
#----------------------------------------
bet = 1/8
gam = 7/4

Tplot = np.linspace(np.min(T), np.max(T), 1024)
SPIN_model = abs(Tplot-Tcrit)**(bet) * 1.09
MS_model = abs(Tplot-Tcrit)**(-gam) * 0.042

#----------------------------------------
plt.figure(figsize=(7,3))


plt.scatter(abs(T-Tcrit), SPIN, marker='x', label = 'Measurements')
plt.plot(abs(Tplot-Tcrit), SPIN_model, ls='--',c='r',label='Power Law Model')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlabel("|T-Tcrit|")
plt.ylabel("Magnetisation Per Spin")
plt.tight_layout()
plt.legend(loc='best')
plt.xlim(0,Tcrit)


#----------------------------------------

plt.figure(figsize=(7,3))

plt.scatter(abs(T-Tcrit), MS, marker='x', label = 'Measurements')
plt.plot(abs(Tplot-Tcrit), MS_model, ls='--',c='r',label='Power Law Model')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlabel("|T-Tcrit|")
plt.ylabel("Magsup Per Lattice Site")
plt.tight_layout()

plt.legend(loc='best')

#----------------------------------------

plt.figure(figsize=(3.5,3))

plt.scatter(T, SPIN, marker='x', label = 'Measurements')
plt.plot(Tplot, SPIN_model, ls='--',c='r',label='Power Law Model')

plt.grid()
plt.xlabel("Temp, T")
plt.ylabel("Magnetisation Per Spin")
plt.tight_layout()
plt.legend(loc='best')
plt.xlim(2,Tcrit)
plt.ylim(min(SPIN),max(SPIN))


#----------------------------------------

plt.figure(figsize=(3.5,3))

plt.scatter(T, MS, marker='x', label = 'Measurements')
plt.plot(Tplot, MS_model, ls='--',c='r',label='Power Law Model')

plt.grid()
plt.xlabel("Temp, T")
plt.ylabel("Magsup Per Lattice Site")
plt.xlim(2,Tcrit)
plt.ylim(min(MS),max(MS))
plt.tight_layout()

plt.legend(loc='best')

#----------------------------------------
plt.show()
