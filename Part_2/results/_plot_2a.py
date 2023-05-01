import numpy as np
import matplotlib.pylab as plt

import warnings
warnings.filterwarnings("ignore")


#----------------------------------------
n        = 32
nsamples = 10000

DATA = np.loadtxt("./temp_sweep_32.dat" )
Tcrit = 2 / np.log(1+np.sqrt(2))

T = DATA[:,0]

CHAINS = DATA[:,-1].astype('int32')

EN, ENVAR = DATA[:,1], DATA[:,2]
SPIN, SPINVAR = DATA[:,3], DATA[:,4]

nchains = len(np.unique(CHAINS))
ntemps = len(CHAINS) // nchains

SPIN /= n**2
EN   /= n**2

SPINVAR /= n**4
ENVAR /= n**4

T1= T[:ntemps]

i_Tcrit = np.where(T1>=Tcrit)[0][0]

#----------------------------------------
# AVERAGE SPIN AND ENERGY OVER ALL CHAINS

SPIN_AV     = np.zeros(ntemps)
SPINVAR_AV  = np.zeros(ntemps)

EN_AV     = np.zeros(ntemps)
ENVAR_AV  = np.zeros(ntemps)


SPINS_BYCHAIN = np.zeros([ntemps,nchains])
SPINVARS_BYCHAIN = np.zeros([ntemps,nchains])

ENS_BYCHAIN = np.zeros([ntemps,nchains])
ENVARS_BYCHAIN = np.zeros([ntemps,nchains])

#===========
for chain in np.unique(CHAINS):
    i=np.where(CHAINS==chain)[0]
    SPINS_BYCHAIN[:,chain]      = SPIN[i]
    SPINVARS_BYCHAIN[:,chain]   = SPINVAR[i]

    ENS_BYCHAIN[:,chain]      = EN[i]
    ENVARS_BYCHAIN[:,chain]   = ENVAR[i]

#===========
for i in range(ntemps):
    if np.min(SPINVARS_BYCHAIN[i,:])!=0:
        weights = SPINVARS_BYCHAIN[i,:]**-2
    else:
        weights = np.ones(nchains)

    SPIN_AV[i]     = np.sum(abs(SPINS_BYCHAIN[i,:])*weights) / np.sum(weights)
    SPINVAR_AV[i]  = np.sum(SPINVARS_BYCHAIN[i,:]*weights) / np.sum(weights)

    #===========

    if np.min(ENVARS_BYCHAIN[i,:])!=0:
        weights = ENVARS_BYCHAIN[i,:]**-2
    else:
        weights = np.ones(nchains)

    EN_AV[i]     = np.sum(ENS_BYCHAIN[i,:]*weights) / np.sum(weights)
    ENVAR_AV[i]  = np.sum(ENVARS_BYCHAIN[i,:]*weights) / np.sum(weights)

#----------------------------------------
# GET HEAT CAPACITY AND MAGNETIC SUSCEPTIBILITY
CP = ENVAR_AV/T1**2
CP_VAR = (np.sum(ENVARS_BYCHAIN**2,axis=1) / nchains -ENVAR_AV**2) /T1**2

MAGSUP = SPINVAR_AV/T1
MAGSUP_VAR = (np.sum(SPINVARS_BYCHAIN**2,axis=1) / nchains -SPINVAR_AV**2) /T1

#----------------------------------------
# PLOT CONFIGURATION
figsize= (7,3)
xlabel = "Temp, $k_B T / J$"
ylabel = "Spin per Lattice Site, $\mathcal{M}$"

Tplot = np.linspace(np.min(T),np.max(T),128)
M_model = np.where(Tplot<Tcrit,(1-np.sinh(2/Tplot)**-4)**(1/8),0)

eta = 1
gam = 7/4

#----------------------------------------
# M VS T, ALL CHAINS

plt.figure(figsize=figsize)

plt.scatter( T, SPIN, marker='.')
plt.errorbar(T, SPIN, yerr = np.sqrt(SPINVAR), fmt='none',alpha=0.1)
plt.errorbar(T, SPIN, yerr = np.sqrt(SPINVAR/nsamples), fmt='none')

plt.plot(Tplot,M_model,c='k',label='Model')
plt.plot(Tplot,-M_model,c='k')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.tight_layout()

plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$",alpha=0.5)
#----------------------------------------
# M VS T, CHAIN WALKS

plt.figure(figsize=figsize)
for chain in np.unique(CHAINS):
    i=np.where(CHAINS==chain)[0]
    plt.errorbar(T[i],SPIN[i],yerr = np.sqrt(SPINVAR[i]),lw=0.5,label="Chain %i" %chain)
    plt.scatter(T[i],SPIN[i])

plt.plot(Tplot,M_model,c='k',label='Model')
plt.plot(Tplot,-M_model,c='k')

plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$",alpha=0.5)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.tight_layout()
        
#----------------------------------------
# M VS T, AVERAGED

plt.figure(figsize=figsize)
plt.scatter( T1, SPIN_AV, marker='.')
plt.errorbar(T1, SPIN_AV, yerr = np.sqrt(SPINVAR_AV/nsamples/nchains), fmt='none', capsize=5, label='Average Spin for Chain')
plt.errorbar(T1, SPIN_AV, yerr = np.sqrt(SPINVAR_AV), fmt=':',alpha=.5, label='Spread Over All Chains')
plt.plot(Tplot,M_model,c='k',label='Model')
plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$",alpha=0.5)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc="best")
plt.tight_layout()

#----------------------------------------
# E VS T, AVERAGED

plt.figure(figsize=figsize)
plt.scatter( T1, EN_AV, marker='.')
plt.errorbar(T1, EN_AV, yerr = np.sqrt(ENVAR_AV/nsamples/nchains), fmt='none', capsize=5, label='Average Spin for Chain')
plt.errorbar(T1, EN_AV, yerr = np.sqrt(ENVAR_AV), fmt=':',alpha=.5, label='Spread Over All Chains')
plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$",alpha=0.5)
plt.axhline(0,ls="-",c='k',lw=2)
plt.xlabel(xlabel)
plt.ylabel("Energy Per Lattice Site")
plt.legend(loc="best")
plt.tight_layout()

#----------------------------------------
# HEAT CAPACITY

plt.figure(figsize=figsize)
plt.scatter( T1, CP, marker='.')
plt.errorbar(T1, CP, yerr = np.sqrt(CP_VAR/nchains), fmt='none', capsize=5, label='Average Heat Capacity')
plt.errorbar(T1, CP, yerr = np.sqrt(CP_VAR), fmt=':',alpha=.5, label='Spread Over All Chains')
plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$",alpha=0.5)
plt.xlabel(xlabel)
plt.axhline(0,c='k',lw=0.5)
plt.xlim(min(T1),max(T1))
plt.ylabel("Heat Capacity Per Lattice Site")
plt.legend(loc="best")
plt.tight_layout()

#----------------------------------------
# MAGNETIC SUSCEPTIBILITY

plt.figure(figsize=figsize)
plt.scatter( T1, MAGSUP, marker='.')
plt.errorbar(T1, MAGSUP, yerr = np.sqrt(MAGSUP_VAR/nchains), fmt='none', capsize=5, label='Average Magnetic Susceptibility')
plt.errorbar(T1, MAGSUP, yerr = np.sqrt(MAGSUP_VAR), fmt=':',alpha=.5, label='Spread Over All Chains')


plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$",alpha=0.5)
plt.xlabel(xlabel)
plt.ylabel("$\chi(T)$, Magsup Per Lattice Site")
plt.axhline(0,c='k',lw=0.5)
plt.xlim(min(T1),max(T1))
plt.legend(loc="best")
plt.tight_layout()

plt.show()


