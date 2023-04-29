import numpy as np
import matplotlib.pylab as plt

import warnings
warnings.filterwarnings("ignore")

#----------------------------------------
figsize= (7,3)
xlabel = "Temp, $k_B T / J$"
ylabel = "Spin per Lattice Site, $\mathcal{M}$"

Tcrit = 2 / np.log(1+np.sqrt(2))
Tplot = np.linspace(0, 5,128)
M_model = np.where(Tplot<Tcrit,(1-np.sinh(2/Tplot)**-4)**(1/8),0)

nsamples = 10000
#----------------------------------------
plt.figure(figsize=figsize)
for n in [8,16,32,64]:


    DATA = np.loadtxt("./temp_sweep_%i.dat" %n)

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

    plt.scatter( T1, SPIN_AV, marker='.')
    plt.errorbar(T1, SPIN_AV, yerr = np.sqrt(SPINVAR_AV/nsamples/nchains), fmt='-', capsize=5, label='L=%i' %n)

plt.grid(True)
plt.plot(Tplot,M_model,c='k',label='Model')
plt.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend(loc="best")
plt.tight_layout()

#----------------------------------------
plt.show()


