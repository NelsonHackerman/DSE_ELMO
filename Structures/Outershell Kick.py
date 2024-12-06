import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve, nsolve
import sympy
global E,v,rho,Sig_tu,Sig_ty,g,SF,m,l,d,r,l_ax_com,l_ax_ten,l_lat,m_bend,l_eq_ten,l_eq_com,f_ax,f_lat
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,m,l,d,r,l_ax_com,l_ax_ten,l_lat,m_bend,l_eq_ten,l_eq_com,f_ax,f_lat

GREEN = "\033[92m"
RED = "\033[91m"
PURPLE = "\033[95m"
RESET = "\033[0m"
#Kickstage option, I'm running the program twice, first for the top half with the orbiter, then for the actual kickstage.

mtot=4317+2665
m=2665



def Shell(m,l,d,t,n):
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
    t3=l_eq_ten*SF/(np.pi*d*Sig_ty)
    print('Estimated thickness tension mode: ',t3*1000,' mm')

    #Using t from the input, as buckling is the most critical and has to be done manually as k has to be read off a graph
    #Sizing for compressive buckling strength with bending moment included
    A=np.pi*d*t
    Req_buck=l_eq_com/A 
    #n=200 #number of stringers
    b=np.pi*d/n #distance between stringers
    ts=0.001 #thickness of 1 stringer
    ls=0.02 # total length of 1 stringer
    val1=b**2/r/t*np.sqrt(1-v**2) #A value to read off k from the graph in SMAD
    k=input(f'Read off graph, r/t: {GREEN}{r/t}{RESET}, x axis: {GREEN}{val1}{RESET}: k=')
    k=float(k)
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
    print(f"{PURPLE}Mass of load-bearing structure: {struc_m} kg{RESET}")

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
    print('Total MOI: ',I_total,'m4')
    print('Required MOI ',I_req,'m4')
    f_lat_mod=0.56*np.sqrt(E*I_total/(m*l**3))
    f_ax_mod=0.25*np.sqrt(A*E/(m*l))
    v_struc=np.pi*r**2*l
    print('Axial natural frequency: ',f_ax_mod,'Hz')
    print('Lateral natural frequency: ',f_lat_mod,'Hz')
    print('Total volume',v_struc,'m3')
    return struc_m,v_struc
m_u_stage,v_u_stage=Shell(m,3,d,0.001,190) #for orbiter stage
m_l_stage,v_l_stage=Shell(mtot,3,d,0.001,185) #for the lower stage
total_struc_mass=m_u_stage+m_l_stage
total_struc_vol=v_u_stage+v_l_stage
print(f"{PURPLE}Total structure mass: {total_struc_mass} kg{RESET}")
print(f"{PURPLE}Total outer volume: {total_struc_vol} m^3{RESET}")