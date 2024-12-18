import numpy as np

from sympy import symbols, Eq, solve, nsolve
import sympy
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,m,l,d,r,l_ax_com,l_ax_ten,l_lat,m_bend,l_eq_ten,l_eq_com,f_ax,f_lat
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

#OPTION 2: shell with longitudinal stiffeners
#USE OPTION B OR D ON CONSTANTS FILE



#sizing for axial natural frequency
A=f_ax**2*m*l/(E*0.25**2) #Calculating minimum area
t1=A/(np.pi*d) #Calculating minimum thickness from area
print('Estimated thickness ax vibration: ',t1*1000,' mm')



#Sizing for longitudinal natural frequency
I_req=f_lat**2*m*l**3/(0.56**2*E) #minimum MOI

print('Required MOI ',I_req)
t2=I_req/(np.pi*r**3) #minimum thickness for MOI, Im ignoring the stiffeners rn because the skin can withstand this alone
print('Estimated thickness lat vibration(only skin): ',t2*1000,' mm')


#Sizing for tensile strength with bending moment included
t3=l_eq_ten*m/(np.pi*d*Sig_ty)
print('Estimated thickness tension mode: ',t3*1000,' mm')

#thickness t1 is the most limiting right now so I will be using that
t=max(t1,t2,t3)

t=0.012
#Sizing for compressive buckling strength with bending moment included
#t=0.001 #Manually changing the thickness here to meet buckling req
A=np.pi*d*t

Req_buck=l_eq_com*m/A 
n=220 #number of stringers
b=np.pi*d/n #distance between stringers
ts=0.008 #thickness of 1 stringer
ls=0.02 # total length of 1 stringer
val1=b**2/r/t*np.sqrt(1-v**2) #A value to read off k from the graph in SMAD
print(f'Read off graph, r/t: {GREEN}{r/t}{RESET}, x axis: {GREEN}{val1}{RESET}')
k=4
Mod_buck=k*np.pi**2*E/(12*(1-v**2))*(t/b)**2
print('Required buckling strength: ',Req_buck,' Pa')
print('Model buckling strength: ',Mod_buck,' Pa')
print('Distance between stringers: ',b*100,' cm')
print('Buckling ratio: ',Mod_buck/Req_buck)
if Mod_buck>Req_buck:
    print(f"{GREEN}Model can withstand buckling{RESET}")
else:
    print(f"{RED}Buckles{RESET}")

struc_m=t*d*np.pi*l*rho + ts*ls*l*n*rho #Calculating the total structural mass
print('Mass of load bearing structure: ',struc_m,' kg')
struc_m=1/0.6347*struc_m
print('Mass of TOTAL structure: ',struc_m,' kg')

#Stringer and skin total MOI calculation
I_skin=np.pi*r**3*t #Skin MOI
angle=2*np.pi/n #angle between stringers
lm=[]
lm.extend([0,0,r,r])
for i in range (1,int(n/4)):
    d=float(np.sin(i*angle)*r)
    lm.extend([d,d,d,d]) #this becomes a list of distances of each stringer from the center axis
I_str=np.sum(np.square(lm)*ls*ts) #Stringer MOI, Only parallel axis, inherent MOI ignored
I_total=I_str+I_skin

III_str=1/12*(ts**3*ls/2)+1/12*(ts*(ls/2)**3)*n
print('Total MOI: ',I_total,'m4')
print('Required MOI ',I_req,'m4')
f_lat_mod=0.56*np.sqrt(E*I_total/(m*l**3))
f_ax_mod=0.25*np.sqrt((A+ts*ls*n)*E/(m*l))
print('Axial natural frequency: ',f_ax_mod,'Hz')
print('Lateral natural frequency: ',f_lat_mod,'Hz')
print('Total volume',np.pi*r**2*l,'m3')
print(m)
print('Fraction of MOI: ',III_str/I_total)