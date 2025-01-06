import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve, nsolve
import sympy
global E,v,rho,Sig_tu,Sig_ty,g,SF,l_ax_com,l_ax_ten,l_lat,f_ax,f_lat
from Constants import E,v,rho,Sig_tu,Sig_ty,g,SF,morb,mtot,l_ax_com,l_ax_ten,l_lat,m_bend_u,m_bend_l,l_eq_ten_u,l_eq_ten_l,l_eq_com_u,l_eq_com_l,f_ax,f_lat, ru,rl,lu,ll,dl,du
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

#option 1: monocoque- no ring or longitudinal stiffeners


def Shell(m,l,d,t,l_eq_ten,l_eq_com,r):
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


    #Sizing for tensile strength 
    t3=l_eq_ten*m/(np.pi*d*Sig_ty)
    print('Estimated thickness tension mode: ',t3*1000,' mm')

    #thickness t1 is the most limiting right now so I will be using that
    #t=max(t1,t2,t3)
    

    #Sizing for compressive strength


    A=np.pi*d*t

    Req_buck=l_eq_com*m/A #Minimum required buckling strength
    # phi=1/16*np.sqrt(r/t)
    # gamma=1-0.901*(1-np.exp(-phi))
    Mod_buck=.6*E*t/r*(1-0.901*(1-np.exp(-1/16*np.sqrt(r/t))))
    print('Required buckling strength: ',Req_buck,' Pa')
    print('Model buckling strength: ',Mod_buck,' Pa')
    print('Buckling ratio: ',Mod_buck/Req_buck)
    if Mod_buck>Req_buck:
        print(f"{GREEN}Model can withstand buckling{RESET}")
    else:
        print(f"{RED}Buckles{RESET}")
    # x=np.linspace(0.0001,0.006,100)
    # y=0.6*E*x/r*(1-0.901*(1-np.exp(-1/16*np.sqrt(r/x))))
    # for i in range (len(x)):
    #         print(x[i],y[i]/Req_buck)
    struc_m=t*d*np.pi*l*rho
    print('Mass of load bearing structure: ',struc_m,' kg')
    v_struc=np.pi*r**2*l
    return struc_m,v_struc
m_u_stage,v_u_stage=Shell(morb,lu,du,0.003,l_eq_ten_u,l_eq_com_u,ru) #for orbiter stage
m_l_stage,v_l_stage=Shell(mtot,ll,dl,0.00,l_eq_ten_l,l_eq_com_l,rl) #for the lower stage