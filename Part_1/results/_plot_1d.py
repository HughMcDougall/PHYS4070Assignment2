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
DATA = np.loadtxt("pt_1d.dat")
N           = (DATA.shape[1]-1)//4
N_asteroids = N-2

T = DATA[:,0]
Xmoon1, Ymoon1 = DATA[:,1], DATA[:,2]
Xmoon2, Ymoon2 = DATA[:,3], DATA[:,4]


#----------------------------
#Plot Trajectories

plt.figure(figsize = (5,5))

plt.plot(Xmoon1, Ymoon1, label = "Moon 1", ls=':')
plt.plot(Xmoon2, Ymoon2, label = "Moon 2", ls=':')

Xstarts = np.zeros(N_asteroids)
Ystarts = np.zeros(N_asteroids)

R = [None]*N_asteroids

for i in range(2,N_asteroids+2):
    X = DATA[:, 2*i + 1]
    Y = DATA[:, 2*i + 2]
    

    Xstarts[i-2], Ystarts[i-2] = X[0], Y[0]
    R[i-2] = (X**2+Y**2)**0.5
    
    if i==2:
        plt.plot(X, Y, c='k', lw=0.1, label = 'Asteroid paths')
    else:
        plt.plot(X, Y, c='k', lw=0.1)
plt.scatter(Xstarts,Ystarts,c='k', label = "Asteroid start locations")

plt.plot(np.cos(theta_grid)*R_planet,np.sin(theta_grid)*R_planet)

plt.legend(loc='best')
plt.axis('equal')
plt.xlim(-25,25)
plt.ylim(-25,25)
#----------------------------
plt.figure(figsize=(5,3))
for i in range(N_asteroids):
    if i ==1:
        plt.plot(T,R[i],lw=0.1,c='k',label="Asteroids")
    else:
        plt.plot(T,R[i],lw=0.1,c='k')

plt.title("Asteroid Orbital Radius vs Time")
plt.ylim(0,19*2)
plt.axhline(19,c='k',ls='--', label='Initial')
plt.xlabel("Time")
plt.ylabel("Distance to Planet")
plt.legend(loc='upper right')
plt.tight_layout()
#----------------------------

plt.show()
