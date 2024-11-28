import numpy as np
import matplotlib.pyplot as plt

########################### user input #################################################
# variables
h_orb = 100e+3
inc = np.pi/2
duration= 3*(1.37*24*3600)
dt = 10
phase_tethys = np.pi/4
########################################################################################

# constants
G = 6.67e-11
M_saturn = 5.6832e+26
M_enc = 1.08e+20
a_enc = 238e+6
r_enc = 257e+3
a_tethys = 294.66e+6
T_enc = 2*np.pi*np.sqrt(a_enc**3/G/M_saturn)
M_tethys = 6.18e+20
T_tethys = 2*np.pi*np.sqrt(a_tethys**3/G/M_saturn)
V_enc = np.sqrt(G*M_saturn/a_enc)


def dist(r1,r2):
    # calculates distance between two position vectors r1 and r2
    return np.sqrt((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2 )


def Fg(M,r1,r2):
    # r1 is the big mass location
    # r2 is the small mass location
    return G*M/dist(r1,r2)**3 * (r1-r2)

def pos_enc(t):
    # position of enceladus as a function of time
    x = a_enc*np.cos(2*np.pi/T_enc*t)
    y = a_enc*np.sin(2*np.pi/T_enc*t)
    z = 0
    return np.array([x,y,z])

def pos_tet(t):
    # position of tethys as a function of time
    x = a_tethys*np.cos(2*np.pi/T_tethys*t + phase_tethys)
    y = a_tethys*np.sin(2*np.pi/T_tethys*t + phase_tethys)
    z = 0
    return np.array([x,y,z])



def V_orb(h):
    # orbital velocity of ELMO around enceladus for a circular orbit
    return np.sqrt(G*M_enc/(r_enc+h))

def V_start(h,inc):
    # starting velocity vector of ELMO in the simulation
    Vmag_rel = V_orb(h)
    Vy_rel = np.cos(inc)*Vmag_rel
    Vz_rel = np.sin(inc)*Vmag_rel
    return np.array([0,Vy_rel + V_enc,Vz_rel])

t = 0


nsteps = int(duration/dt)
pos_elmo = np.array([a_enc + r_enc + h_orb,0,0])
vel_elmo = V_start(h_orb, inc)
pos_saturn = np.array([0,0,0])
poslist_elmo = np.zeros((3,nsteps))
vellist_elmo = np.zeros((3,nsteps))
poslist_enc = np.zeros((3,nsteps))
relposlist = np.zeros((3,nsteps))
distlist = np.zeros(nsteps)
timelist = np.arange(0,nsteps*dt,dt)
enceladus_surface = [[0,nsteps*dt/3600],[0, 0]]
for i in range(nsteps):
    
    pos_enceladus = pos_enc(t)
    pos_tethys = pos_tet(t)
    poslist_enc[:,i] = pos_enceladus
    poslist_elmo[:,i] = pos_elmo
    vellist_elmo[:,i] = vel_elmo
    relposlist[:,i] = pos_elmo - pos_enceladus
    distlist[i] = dist(pos_elmo,pos_enceladus)
    
    # the actual physics
    acc = Fg(M_saturn,pos_saturn,pos_elmo) + Fg(M_enc,pos_enceladus,pos_elmo) #+ Fg(M_tethys,pos_tethys, pos_elmo)
    
    vel_elmo += acc*dt
    pos_elmo += vel_elmo*dt
    
    t += dt


#plt.plot(poslist_enc[0,:],poslist_enc[1,:])
#plt.plot(poslist_elmo[0,:],poslist_elmo[1,:])
#plt.plot(relposlist[0,:],relposlist[1,:])
plt.plot(timelist/3600,(distlist-r_enc)/1000)
plt.plot(enceladus_surface[0],enceladus_surface[1])
plt.xlabel("Time [hr]")
plt.ylabel("Altitude [km]")

plt.grid()
#plt.legend(['Enceladus','ELMO'])

