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

r_soi    = r_moon / (1+np.sqrt(10))

theta_grid = np.linspace(0,2*np.pi,512)

#----------------------------
#LOAD DATA
DATA_A = np.loadtxt("pt_1a.dat")
DATA_B = np.loadtxt("pt_1b.dat")

T_a = DATA_A[:,0]
T_b = DATA_B[:,0]

Xproj_a, Yproj_a = DATA_A[:,1], DATA_A[:,2]
Xproj_b, Yproj_b = DATA_B[:,1], DATA_B[:,2]

Xmoon_a, Ymoon_a = DATA_A[:,3],   DATA_A[:,4]
Xmoon_b, Ymoon_b = DATA_B[:,3],   DATA_B[:,4]

dR_a = ((Xproj_a-Xmoon_a)**2 + (Yproj_a-Ymoon_a)**2)**0.5
dR_b = ((Xproj_b-Xmoon_b)**2 + (Yproj_b-Ymoon_b)**2)**0.5


#Itteration no's of max altitude and closest approach
i_a = find_nearest(dR_a, np.min(dR_a))
i_b = find_nearest(dR_b, R_moon)

i_soi= np.where(dR_b<r_soi)[0][0]

#----------------------------
#Plot Trajectories

plt.figure(figsize = (5,5))

plt.plot(Xproj_a[:i_b], Yproj_a[:i_b], label = "Trajectory, Planet Only", ls='--', c='k', lw=0.5)
plt.plot(Xproj_b[:i_b], Yproj_b[:i_b], label = "Trajectory, Planet + Moon")
plt.plot(Xmoon_b[:i_b], Ymoon_b[:i_b], label = "Moon Path", ls=':')
plt.plot(Xmoon_a, Ymoon_a,                                  ls=':', c='k', alpha=0.5)

plt.plot(np.cos(theta_grid)*R_planet,np.sin(theta_grid)*R_planet)

plt.scatter(Xproj_a[i_b], Yproj_a[i_b], c='k', alpha=0.5)
plt.scatter(Xproj_b[i_b], Yproj_b[i_b])
plt.scatter(Xmoon_b[i_b], Ymoon_b[i_b])
plt.plot(np.cos(theta_grid)*R_moon+Xmoon_b[i_b],np.sin(theta_grid)*R_moon+Ymoon_b[i_b])

plt.scatter(Xproj_b[i_soi], Yproj_b[i_soi], c='k', marker='x', label = 'Projectile Enters Moon SOI')


plt.scatter(Xmoon_a[i_a], Ymoon_a[i_a], c='k', alpha=0.5)
plt.plot(np.cos(theta_grid)*R_moon+Xmoon_a[i_a],np.sin(theta_grid)*R_moon+Ymoon_b[i_a], c='k', alpha=0.5)

plt.legend(loc='best')
plt.axis('equal')


#----------------------------
#Plot path convergence


fig, ax = plt.subplots(2,1, figsize=(5,4), sharex=True)

ax[0].set_title("Projectile / Moon Distance")
ax[0].plot(T_a[:i_b], dR_a[:i_b], ls='--', c='k', lw=0.5)
ax[0].plot(T_b[:i_b], dR_b[:i_b])
ax[0].axhline(R_moon,ls='--')
ax[0].axhline(r_soi,ls=':',label="Moon Sphere of Influence", c='k')
ax[0].scatter(T_b[i_soi], dR_b[i_soi], c='k', marker='x', label = 'Projectile Enters Moon SOI')
ax[0].legend(loc='best')

ax[1].set_title("Projectile Altitude")
ax[1].plot(T_a[:i_b], ((Xproj_a**2+Yproj_a**2)**0.5)[:i_b], label = "Trajectory, Planet Only", ls='--', c='k', lw=0.5)
ax[1].plot(T_b[:i_b], ((Xproj_b**2+Yproj_b**2)**0.5)[:i_b], label = "Trajectory, Planet + Moon")

ax[1].axhline(r_moon,ls='--', c='k',lw=0.5)
ax[1].axhline(r_moon-R_moon,ls='--', c='k',lw=0.5)
ax[1].legend(loc='best')

ax[1].set_xlabel("Time")

fig.tight_layout()

print("Time of impact for new trajectory:\t t = %.2f" %T_a[i_b])

#----------------------------

plt.show()
