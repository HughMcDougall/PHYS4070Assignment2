import numpy as np
import matplotlib.pylab as plt

import warnings
warnings.filterwarnings("ignore")


#----------------------------------------

nsamples = 10000
Tcrit = 2 / np.log(1+np.sqrt(2))

#----------------------------------------
# PLOT CONFIGURATION
figsize= (7,3)
xlabel = "Temp, $k_B T / J$"
ylabel = "Spin per Lattice Site, $\mathcal{M}$"

Tplot = np.linspace(0,5,128)
M_model = np.where(Tplot<Tcrit,(1-np.sinh(2/Tplot)**-4)**(1/8),0)

nu = 1.0
gam = 7.0/4.0

ZING=[]
#----------------------------------------
MFIG, MAX = plt.subplots(1,1, figsize=figsize)
EFIG, EAX = plt.subplots(1,1, figsize=figsize)
CPFIG, CPAX = plt.subplots(1,1, figsize=figsize)
MSFIG, MSAX = plt.subplots(1,1, figsize=figsize)
CPFIG_scale, CPAX_scale = plt.subplots(1,1, figsize=figsize)
MSFIG_scale, MSAX_scale = plt.subplots(1,1, figsize=figsize)

for n in [8,16,32,64]:
    print(n)

    DATA = np.loadtxt("./temp_sweep_%i.dat" %n)

    T = DATA[:,0]

    CHAINS = DATA[:,-1].astype('int32')

    EN, ENVAR = DATA[:,1], DATA[:,2]
    SPIN, SPINVAR = DATA[:,3], DATA[:,4]

    nchains = len(np.unique(CHAINS))
    ntemps = len(CHAINS) // nchains

    if True:
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


    if True:
        CP*=n**2
        CP_VAR*=n**4
        
        MAGSUP*=n**2
        MAGSUP_VAR*=n**4

    fac = np.log(n)
    CP_scale=CP/fac
    CP_VAR_scale=CP_VAR/fac**2

    fac = n**(gam/nu)
    MAGSUP_scale=MAGSUP/fac
    MAGSUP_VAR_scale=MAGSUP_VAR/fac**2

    Tscale = (T1/Tcrit-1)*n**(1/nu)

    #----------------------------------------
    #PLOTS
    #----------------------------------------
    # M VS T, AVERAGED
    MAX.scatter( T1, SPIN_AV, marker='.')
    MAX.errorbar(T1, SPIN_AV, yerr = np.sqrt(SPINVAR_AV/nsamples/nchains), fmt='-', capsize=5, label='L = %i' %n)

    #----------------------------------------
    # E VS T, AVERAGED
    EAX.scatter( T1, EN_AV, marker='.')
    EAX.errorbar(T1, EN_AV, yerr = np.sqrt(ENVAR_AV/nsamples/nchains), fmt='-', capsize=5, label='L = %i' %n)
    EAX.errorbar(T1, EN_AV, yerr = np.sqrt(ENVAR_AV), fmt='-', alpha=0.5)

    #----------------------------------------
    # HEAT CAPACITY

    CPAX.scatter( T1, CP, marker='.')
    CPAX.errorbar(T1, CP, yerr = np.sqrt(CP_VAR/nchains), fmt='-', capsize=5, label='L = %i' %n)

    #----------------------------------------
    # MAGNETIC SUSCEPTIBILITY

    MSAX.scatter( T1, MAGSUP, marker='.')
    MSAX.errorbar(T1, MAGSUP, yerr = np.sqrt(MAGSUP_VAR/nchains), fmt='-', capsize=5, label='L = %i' %n)

        #----------------------------------------
    # HEAT CAPACITY

    CPAX_scale.scatter( Tscale, CP_scale, marker='.')
    CPAX_scale.errorbar(Tscale, CP_scale, yerr = np.sqrt(CP_VAR_scale/nchains), fmt='-', capsize=5, label='L = %i' %n)

    #----------------------------------------
    # MAGNETIC SUSCEPTIBILITY

    MSAX_scale.scatter( Tscale, MAGSUP_scale, marker='.')
    MSAX_scale.errorbar(Tscale, MAGSUP_scale, yerr = np.sqrt(MAGSUP_VAR_scale/nchains), fmt='-', capsize=5, label='L = %i' %n)


#----------------------------------------
# M VS T, AVERAGED
MAX.plot(Tplot,M_model,c='k',label='Model')
MAX.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
MAX.set_xlabel(xlabel)
MAX.set_ylabel(ylabel)
MAX.legend(loc="best")
MFIG.tight_layout()

#----------------------------------------
# E VS T, AVERAGED
EAX.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
EAX.axhline(0,ls="-",c='k',lw=2)
EAX.set_xlabel(xlabel)
EAX.set_ylabel("Energy Per Lattice Site")
EAX.legend(loc="best")
EFIG.tight_layout()
#----------------------------------------
# HEAT CAPACITY

CPAX.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
CPAX.set_xlabel(xlabel)
CPAX.set_ylabel("$c_p$ Per Site")
CPAX.legend(loc="best")
#CPAX.set_xlim(0,max(Tspace))
#CPAX.set_yscale("log")
#CPAX.set_xscale("log")
CPAX.grid()
CPFIG.tight_layout()

#----------------------------------------
# MAGNETIC SUSCEPTIBILITY

MSAX.axvline(Tcrit,ls='--',c='k',label="$T_{crit}$")
MSAX.set_xlabel(xlabel)
MSAX.set_ylabel("$Mag Sup Per Site$")
MSAX.legend(loc="best")
#MSAX.set_xlim(0,max(Tspace))
#MSAX.set_yscale("log")
#MSAX.set_xscale("log")
MSAX.grid()
MSFIG.tight_layout()

#----------------------------------------
# HEAT CAPACITY_SCALE

CPAX_scale.axvline(0,ls='--',c='k',label="$T_{crit}$")
CPAX_scale.set_xlabel("$(T/T_{crit}-1)L^{1/\\nu }$")
CPAX_scale.set_ylabel("$c_p / log_(L)$")
CPAX_scale.legend(loc="best")
#CPAX.set_xlim(0,max(Tspace))
#CPAX.set_yscale("log")
#CPAX.set_xscale("log")
CPAX_scale.grid()
CPFIG_scale.tight_layout()

#----------------------------------------
# MAGNETIC SUSCEPTIBILITY_SCALE

MSAX_scale.axvline(0,ls='--',c='k',label="$T_{crit}$",lw=0.2)
MSAX_scale.set_xlabel("$(T/T_{crit}-1)L^{1/\\nu}$")
MSAX_scale.set_ylabel("$M / L^{\gamma / \\nu}$")
MSAX_scale.legend(loc="best")
#MSAX.set_xlim(0,max(Tspace))
MSAX_scale.set_yscale("log")
#MSAX.set_xscale("log")
MSAX_scale.grid()
MSFIG_scale.tight_layout()

#----------------------------------------
plt.show()
