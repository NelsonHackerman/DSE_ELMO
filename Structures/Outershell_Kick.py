import numpy as np

from sympy import symbols, Eq, solve, nsolve
import sympy
global E,v,rho,Sig_tu,Sig_ty,g,SF,l_ax_com,l_ax_ten,l_lat,f_ax,f_lat,Sig_c
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,morb,mtot,l_ax_com,l_ax_ten,l_lat,m_bend_u,m_bend_l,l_eq_ten_u,l_eq_ten_l,l_eq_com_u,l_eq_com_l,f_ax,f_lat, ru,rl,lu,ll,dl,du,Sig_c

GREEN = "\033[92m"
RED = "\033[91m"
PURPLE = "\033[95m"
RESET = "\033[0m"
#Kickstage option, I'm running the program twice, first for the top half with the orbiter, then for the actual kickstage.




def Shell(m,l,d,t,n,l_eq_ten,l_eq_com,r):
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
    #Using t from the input, as buckling is the most critical and has to be done manually as k has to be read off a graph
    #Sizing for compressive buckling strength with bending moment included
    A=np.pi*d*t
    realA=np.pi*((r+t/2)**2-(r-t/2)**2) #checking accuraacy of thin wall assump
    Req_buck=l_eq_com*m/A 
    realReq_buck=l_eq_com*m/realA
    buckratio=(realReq_buck-Req_buck)/Req_buck*100 #thin wall assump check
    
    b=np.pi*d/n #distance between stringers
    ts=0.001 #thickness of 1 stringer
    ls=0.03 # total length of 1 stringer
    A2=np.pi*d*t+ls*ts*n
    if l_eq_com*m/A2>Sig_c:
        print(f"{RED}Compression{RESET}")
        print(l_eq_com*m/A2/Sig_c)
    else:
        print(f"{GREEN}No compression{RESET}")
        print(l_eq_com*m/A2/Sig_c)
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
    m_shell=t*d*np.pi*l*rho
    print('Mass of shell ONLY: ',m_shell)

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
    f_ax_mod=0.25*np.sqrt((A+ls*ts*n)*E/(m*l))
    v_struc=np.pi*r**2*l
    print('Axial natural frequency: ',f_ax_mod,'Hz')
    print('Lateral natural frequency: ',f_lat_mod,'Hz')
    print('Total volume',v_struc,'m3')
    print('Percentage difference assumptions',(realA-A)/A*100)
    print('Check assumptions',r/t,l/r)
    print('Assumption ratio',buckratio)
    print('End of calc')
    return struc_m,v_struc
m_u_stage,v_u_stage=Shell(morb,lu,du,0.001,104,l_eq_ten_u,l_eq_com_u,ru) #for orbiter stage
m_l_stage,v_l_stage=Shell(mtot,ll,dl,0.001,208,l_eq_ten_l,l_eq_com_l,rl) #for the lower stage
total_struc_mass=m_u_stage+m_l_stage
total_struc_vol=v_u_stage+v_l_stage
print(f"{PURPLE}Total LOAD BEARING mass: {total_struc_mass} kg{RESET}")
total_struc_mass=1/0.6347*total_struc_mass #1.1 for fasteners


print(f"{PURPLE}Total structure mass: {total_struc_mass} kg{RESET}")
print(f"{PURPLE}Total outer volume: {total_struc_vol} m^3{RESET}")
print(m_u_stage)
print(m_l_stage)
m_u_stage=1/0.6347*m_u_stage
m_l_stage=1/0.6347*m_l_stage
print ('upper stage mass ', m_u_stage)
print('lower stage mass ',m_l_stage)

