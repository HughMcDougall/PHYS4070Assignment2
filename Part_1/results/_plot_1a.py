import numpy as np
import matplotlib.pylab as plt
#----------------------------
def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx
#----------------------------
#SYSTEM PROPS
R_planet = 1
R_moon   = 0.25
r_moon   = 19

theta_grid = np.linspace(0,2*np.pi,512)

#----------------------------
#LOAD DATA
DATA = np.loadtxt("pt_1a.dat")

T = DATA[:,0]
Xproj, Yproj = DATA[:,1],DATA[:,2]
Xmoon, Ymoon = DATA[:,3],DATA[:,4]
dR = ((Xproj-Xmoon)**2 + (Yproj-Ymoon)**2)**0.5

#Itteration no's of max altitude and closest approach
i_orb       = find_nearest(Xmoon, np.max(Xmoon))
i_nearest   = find_nearest(dR, np.min(dR))

#----------------------------
#Plot Trajectories

plt.figure(figsize = (5,5))

plt.plot(Xproj, Yproj, label = "Projectile Path")
plt.plot(Xmoon, Ymoon, label = "Moon Path", ls=':')

plt.plot(np.cos(theta_grid)*R_planet,np.sin(theta_grid)*R_planet)

plt.scatter(Xproj[i_orb], Yproj[i_orb], label = "Projectile Position at Max Altitude")
plt.scatter(Xmoon[i_orb], Ymoon[i_orb], label = "Moon Position at Max Altitude")
plt.plot(np.cos(theta_grid)*R_moon+Xmoon[i_orb],np.sin(theta_grid)*R_moon+Ymoon[i_orb])

plt.scatter(Xproj[i_nearest], Yproj[i_nearest], label = "Projectile Position at Closest Approach")
plt.scatter(Xmoon[i_nearest], Ymoon[i_nearest], label = "Moon Position at Closest Approach")
plt.plot(np.cos(theta_grid)*R_moon+Xmoon[i_nearest],np.sin(theta_grid)*R_moon+Ymoon[i_nearest])

plt.legend(loc='best')
plt.axis('equal')

print("Maximum altitude is x=%.2f at time t=%.2f" %(Xproj[i_orb], T[i_orb]))


#----------------------------
#Plot path convergence

fig, ax = plt.subplots(2,1, figsize=(5,4), sharex=True)

ax[0].set_title("Projectile / Moon Distance")
ax[0].plot(T,dR)
ax[0].axhline(R_moon,ls='--')
ax[0].axvline(T[i_orb],ls=':')
ax[0].axvline(T[i_nearest],ls=':')

ax[1].set_title("Projectile Altitude")
ax[1].plot(T,(Xproj**2+Yproj**2**0.5))

ax[1].axhline(r_moon,ls='--', c='k',lw=0.5)
ax[1].axhline(r_moon-R_moon,ls='--', c='k',lw=0.5)

ax[1].axvline(T[i_orb],ls=':',c='k')
ax[1].set_xlabel("Time")


fig.tight_layout()
print("Closest approach is |dr|=%.2f at time t=%.2f" %(dR[i_nearest], T[i_nearest]))

#----------------------------
#Plot energy
fig2, ax2 = plt.subplots(2,1,figsize=(5,3))
VprojX, VprojY = DATA[:,5], DATA[:,6]
VmoonX, VmoonY = DATA[:,7], DATA[:,8]

R_proj = (Xproj**2+Yproj**2)**0.5
R_moon = (Xmoon**2+Ymoon**2)**0.5

Eproj = (VprojX**2+VprojY**2)/2 - 1/R_proj
Emoon = (VmoonX**2+VmoonY**2)/2 - 1/R_moon

ax2[0].plot(T,abs(Eproj/Eproj[0]-1))
ax2[0].set_yscale('log')

ax2[1].plot(T,abs(Emoon/Emoon[0]-1))
ax2[1].set_yscale('log')
ax2[0].set_ylim(1E-9,1E-4)
ax2[1].set_ylim(1E-9,1E-4)

ax2[0].set_xlabel("Time")

ax2[0].set_title("Relative Energy Change")
plt.tight_layout()


plt.show()
