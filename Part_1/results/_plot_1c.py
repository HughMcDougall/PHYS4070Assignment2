import numpy as np
import matplotlib.pylab as plt
from math import atan2
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
DATA = np.loadtxt("pt_1c.dat")

T = DATA[:,0]
Xproj, Yproj = DATA[:,1],DATA[:,2]
Xmoon, Ymoon = DATA[:,3],DATA[:,4]
dR = ((Xproj-Xmoon)**2 + (Yproj-Ymoon)**2)**0.5

#Itteration no's of max altitude and closest approach
aX = Xproj[2:] -2 * Xproj[1:-1] + Xproj[:-2] 
aY = Yproj[2:] -2 * Yproj[1:-1] + Yproj[:-2]
dt = T[1] - T[0]
A  = (aX**2+aY**2)**0.5/dt**2
i_kick       = np.where(A==np.max(A))[0][0] + 2

#----------------------------
#Plot Trajectories

plt.figure(figsize = (5,5))

plt.plot(Xproj, Yproj, label = "Projectile Path")
plt.plot(Xmoon, Ymoon, label = "Moon Path", ls=':')

plt.plot(np.cos(theta_grid)*R_planet,np.sin(theta_grid)*R_planet)

plt.scatter(Xproj[i_kick], Yproj[i_kick], label = "Projectile Position at Kick")
plt.scatter(Xmoon[i_kick], Ymoon[i_kick], label = "Moon Position at Kick")
plt.plot(np.cos(theta_grid)*R_moon+Xmoon[i_kick],np.sin(theta_grid)*R_moon+Ymoon[i_kick])

plt.legend(loc='best')
plt.axis('equal')

V_pre_kick  = np.array([Xproj[i_kick-1] - Xproj[i_kick-2],Yproj[i_kick-1] - Yproj[i_kick-2]]) / dt
V_post_kick = np.array([Xproj[i_kick+2] - Xproj[i_kick+1],Yproj[i_kick+2] - Yproj[i_kick+1]]) / dt
dV = (V_post_kick - V_pre_kick)
dVnorm = np.sqrt(dV[0]**2+dV[1]**2)
angle_pre_kick = atan2(V_pre_kick[1],V_pre_kick[0])
angle_post_kick = atan2(V_post_kick[1],V_post_kick[0])
angle_kick = atan2(dV[1],dV[0])
rel_angle = angle_kick-angle_pre_kick


print("Kick at time t=%.2f. Approx Magnitude = %.2f at angle %.2f° to prograde"    %(T[i_kick], dVnorm, rel_angle * 180/np.pi) )

#----------------------------
#Plot path convergence

fig, ax = plt.subplots(2,1, figsize=(5,4), sharex=True)

ax[0].set_title("Projectile / Moon Distance")
ax[0].plot(T,dR)
ax[0].axhline(R_moon,ls='--')
ax[0].axvline(T[i_kick],ls=':')

ax[1].set_title("Projectile Altitude")
ax[1].plot(T,Xproj)

ax[1].axhline(r_moon,ls='--', c='k',lw=0.5)
ax[1].axhline(r_moon+R_moon * 2,ls=':', c='k',lw=0.5)
ax[1].axhline(r_moon-R_moon * 2,ls=':', c='k',lw=0.5)

ax[1].axvline(T[i_kick],ls=':',c='k')
ax[1].set_xlabel("Time")


fig.tight_layout()


#----------------------------

plt.show()
