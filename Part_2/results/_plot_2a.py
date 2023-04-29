import numpy as np
import matplotlib.pylab as plt

import warnings
warnings.filterwarnings("ignore")


#----------------------------------------
n        = 16
nsamples = 10000

DATA = np.loadtxt("./temp_sweep_%i.dat" %n)
Tcrit = 2 / np.log(1+np.sqrt(2))

T = DATA[:,0]

CHAINS = DATA[:,-1].astype('int32')

EN, ENVAR = DATA[:,1], DATA[:,2]
SPIN, SPINVAR = DATA[:,3], DATA[:,4]

nchains = len(np.unique(CHAINS))
ntemps = len(CHAINS) // nchains

SPIN /= n**2
SPINVAR /= n**4
T1= T[:ntemps]

#----------------------------------------

SPIN_AV     = np.zeros(ntemps)
SPINVAR_AV  = np.zeros(ntemps)

SPINS_BYCHAIN = np.zeros([ntemps,nchains])
SPINVARS_BYCHAIN = np.zeros([ntemps,nchains])
for chain in np.unique(CHAINS):
    i=np.where(CHAINS==chain)[0]
    SPINS_BYCHAIN[:,chain]      = SPIN[i]
    SPINVARS_BYCHAIN[:,chain]   = SPINVAR[i]

for i in range(ntemps):
    if np.min(SPINVARS_BYCHAIN[i,:])!=0:
        weights = SPINVARS_BYCHAIN[i,:]**-2
    else:
        weights = np.ones(nchains)

    SPIN_AV[i]     = np.sum(abs(SPINS_BYCHAIN[i,:])*weights) / np.sum(weights)
    SPINVAR_AV[i]  = np.sum(SPINVARS_BYCHAIN[i,:]*weights) / np.sum(weights)

#----------------------------------------
figsize= (7,3)
xlabel = "Temp, $k_B T / J$"
ylabel = "Spin per Lattice Site, $\mathcal{M}$"

Tplot = np.linspace(np.min(T),np.max(T),128)
M_model = np.where(Tplot<Tcrit,(1-np.sinh(2/Tplot)**-4)**(1/8),0)

#----------------------------------------
plt.figure(figsize=figsize)

plt.scatter( T, SPIN, marker='.')
plt.errorbar(T, SPIN, yerr = np.sqrt(SPINVAR), fmt='none',alpha=0.1)
plt.errorbar(T, SPIN, yerr = np.sqrt(SPINVAR/nsamples), fmt='none')

plt.plot(Tplot,M_model,c='k',label='Model')
plt.plot(Tplot,-M_model,c='k')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.tight_layout()

plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
#----------------------------------------
plt.figure(figsize=figsize)
for chain in np.unique(CHAINS):
    i=np.where(CHAINS==chain)[0]
    plt.errorbar(T[i],SPIN[i],yerr = np.sqrt(SPINVAR[i]),lw=0.5,label="Chain %i" %chain)
    plt.scatter(T[i],SPIN[i])

plt.plot(Tplot,M_model,c='k',label='Model')
plt.plot(Tplot,-M_model,c='k')

plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.tight_layout()
        
#----------------------------------------
plt.figure(figsize=figsize)
plt.scatter( T1, SPIN_AV, marker='.')
plt.errorbar(T1, SPIN_AV, yerr = np.sqrt(SPINVAR_AV/nsamples/nchains), fmt='none', capsize=5, label='Average Spin +- Uncertainty')
plt.errorbar(T1, SPIN_AV, yerr = np.sqrt(SPINVAR_AV), fmt=':',alpha=.5, label='Spread Over All Chains')
plt.plot(Tplot,M_model,c='k',label='Model')
plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc="best")
plt.tight_layout()

#----------------------------------------
plt.show()


