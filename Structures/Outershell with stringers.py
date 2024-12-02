import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve, nsolve
import sympy

#option 1: monocoque- no ring or longitudinal stiffeners
#material: AL 7075
E=71*10**9#N/m2
v=0.33
rho=2.8*10**3#kg/m3
Sig_tu=524*10**6#N/m2 ultimate tensile strength
Sig_ty=448*10**6#N/m2 yield tensile strength
g=9.81
SF=1.6 #table 11-54 SMAD
#falcon 9: axial -2 to 6g, long -2 to 2g + for compress
#arianes: axial -6 to 2.5g, long -1.8 to 1.8g - for compress
#overall axial -6 to 2.5g, long -2 to 2g (I'M USING + FOR TENSION)
#falcon 9: axial >25Hz, lat >10Hz
#arianes: axial >20Hz, lat >6Hz
#overall axial >25Hz, lat >10Hz

m=7355.22 #estimated wet mass from the mass budget
l=7
d=2
r=d/2
l_ax_com=6*g*SF*m
l_ax_ten=2.5*g*SF*m
l_lat=2*g*SF*m
m_bend=l/2*l_lat
f_ax=25
f_lat=10
p=0#pressure for tanks


#sizing for natural frequency

#l=np.sqrt(E/rho)*0.25/f_ax #equation for axial natural frequency
#print('Required length for natural frequency monocoque : ',l)
A=f_ax**2*m*l/(E*0.25**2)

t1=A/(np.pi*d)
print('Estimated thickness ax vibration: ',t1*1000,' mm')



#d=f_lat*l**2/0.56*np.sqrt(8*rho/E)#equation for lateral natural frequency
I=f_lat*m*l**3/(0.56**2*E)
t2=I/(np.pi*r**3)
print('Estimated thickness lat vibration: ',t2*1000,' mm')




#equivalent axial tensile load with bending
l_eq_ten=l_ax_ten+2*m_bend/r

#Sizing for tensile strength 
t3=l_eq_ten*SF/(np.pi*d*Sig_ty)
print('Estimated thickness tension mode: ',t3*1000,' mm')

#thickness t1 is the most limiting right now so I will be using that
t=max(t1,t2,t3)
t=4.4*10**(-3)

#Sizing for compressive strength
l_eq_com=l_ax_com+2*m_bend/r

A=np.pi*d*t

Req_buck=l_eq_com/A #Minimum required buckling strength
# phi=1/16*np.sqrt(r/t)
# gamma=1-0.901*(1-np.exp(-phi))
Mod_buck=.6*E*t/r*(1-0.901*(1-np.exp(-1/16*np.sqrt(r/t))))
print('Required buckling strength: ',Req_buck,' Pa')
print('Model buckling strength: ',Mod_buck,' Pa')
if Mod_buck>Req_buck:
    print('Model can withstand buckling')
else:
    print('Buckles')
# x=np.linspace(0.0001,0.006,100)
# y=0.6*E*x/r*(1-0.901*(1-np.exp(-1/16*np.sqrt(r/x))))
# for i in range (len(x)):
#         print(x[i],y[i]/Req_buck)
struc_m=t*d*np.pi*l*rho
print('Mass of load bearing structure: ',struc_m,' kg')