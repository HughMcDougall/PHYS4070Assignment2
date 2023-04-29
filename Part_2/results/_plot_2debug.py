import numpy as np
import matplotlib.pylab as plt

#----------------------------------------
DATA = np.loadtxt("./temp_sweep_8.dat")
nsamples = 1000
n=8

Tcrit = 2 / np.log(1+np.sqrt(2))

T = DATA[:,0]
CHAINS = DATA[:,-1]

EN, ENVAR = DATA[:,1], DATA[:,2]
SPIN, SPINVAR = DATA[:,3], DATA[:,4]

#----------------------------------------
plt.figure()
plt.scatter(T,SPIN, marker='.')
plt.errorbar(T,SPIN,yerr = np.sqrt(SPINVAR), fmt='none')

Tplot = np.linspace(np.min(T),np.max(T),128)
M_model = np.where(Tplot<Tcrit,(1-np.sinh(2/Tplot)**-4)**(1/8) * n**2,0)
plt.plot(Tplot,M_model,c='k',label='Model')
plt.plot(Tplot,-M_model,c='k')
plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
#----------------------------------------
plt.figure()
for chain in np.unique(CHAINS):
    i=np.where(CHAINS==chain)[0]
    plt.errorbar(T[i],SPIN[i],yerr = np.sqrt(SPINVAR[i]),lw=0.5)
    plt.scatter(T[i],SPIN[i])

plt.plot(Tplot,M_model,c='k',label='Model')
plt.plot(Tplot,-M_model,c='k')
plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
#----------------------------------------
plt.show()


