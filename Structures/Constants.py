#material: AL 7075
# E=71*10**9#N/m2
# v=0.33
# rho=2.8*10**3#kg/m3
# Sig_tu=524*10**6#N/m2 ultimate tensile strength
# Sig_ty=448*10**6#N/m2 yield tensile strength

#t300
#E=140*10**9#N/m2
#v=0.32
#rho=1.636*10**3#kg/m3
#Sig_tu=1820*10**6#N/m2 ultimate tensile strength
#Sig_c=1470*10**6#N/m2 compressive strength
#Sig_ty=1820*10**6#N/m2 yield tensile strength

#msx
E=9.28*10**10#N/m2
v=0.33
rho=1.636*10**3#kg/m3 correct is 1661
Sig_tu=3.46*10**8#N/m2 ultimate tensile strength
Sig_c=3.34*10**8#N/m2 compressive strength
Sig_ty=Sig_tu#N/m2 yield tensile strength



#print(Sig_tu/1.25)
#print(Sig_ty/1.1)#Yeild stress is limiting here
# E=1.150*10**9#N/m2
# v=0.3
# rho=130#kg/m3
# Sig_tu=9.9*10**6#N/m2 ultimate tensile strength
# Sig_ty=9.9*10**6
# E=193*10**9#N/m2
# v=0.3
# rho=8*10**3#kg/m3
# Sig_tu=515*10**6#N/m2 ultimate tensile strength
# Sig_ty=205*10**6
g=9.81
SF=1.25 #table 11-54 SMAD
#falcon 9: axial -2 to 6g, long -2 to 2g + for compress
#arianes: axial -6 to 2.5g, long -1.8 to 1.8g - for compress
#overall axial -6 to 2.5g, long -2 to 2g (I'M USING + FOR TENSION)
#falcon 9: axial >25Hz, lat >10Hz
#arianes: axial >20Hz, lat >6Hz
#overall axial >25Hz, lat >10Hz
l_ax_com=6*g*SF #axial load for compression
l_ax_ten=2.5*g*SF #axial load for tension
l_lat=1.8*g*SF #lateral load
 #axial compressive load with bending
f_ax=20 #axial natural frequency requirement
f_lat=6 #longitudinal natural frequency requirement

option='c'
if option=='a':
    morb=4364
    mtot=10761
    ru=0.75
    rl=1.04
    du=2*ru
    dl=2*rl
    lu=3.35
    ll=1.75
    m_bend_u=lu/2*l_lat #axial bending moment from lateral loads
    l_eq_ten_u=l_ax_ten+2*m_bend_u/ru #axial tension load with bending
    l_eq_com_u=l_ax_com+2*m_bend_u/ru
    m_bend_l=(lu+ll)/2*l_lat#axial bending moment from lateral loads
    l_eq_ten_l=l_ax_ten+2*m_bend_u/rl #axial tension load with bending
    l_eq_com_l=l_ax_com+2*m_bend_u/rl
if option=='b':
    m=17696
    r=1.04
    l=4.55
    d=2*r
    m_bend=l/2*l_lat #axial bending moment from lateral loads
    l_eq_ten=l_ax_ten+2*m_bend/r #axial tension load with bending
    l_eq_com=l_ax_com+2*m_bend/r
if option=='c':
    morb=4308.30
    mtot=6307.75+morb
    ru=1.503/2
    rl=1.944/2
    du=2*ru
    dl=2*rl
    lu=2.55+0.09
    ll=4.25-0.078
    m_bend_u=lu/2*l_lat #axial bending moment from lateral loads
    l_eq_ten_u=l_ax_ten+2*m_bend_u/ru #axial tension load with bending
    l_eq_com_u=l_ax_com+2*m_bend_u/ru
    m_bend_l=(lu+ll)/2*l_lat#axial bending moment from lateral loads
    l_eq_ten_l=l_ax_ten+2*m_bend_l/rl #axial tension load with bending
    l_eq_com_l=l_ax_com+2*m_bend_l/rl
if option=='d':
    m=31047 
    r=1.04
    l=7.54
    d=2*r
    m_bend=l/2*l_lat #axial bending moment from lateral loads
    l_eq_ten=l_ax_ten+2*m_bend/r #axial tension load with bending
    l_eq_com=l_ax_com+2*m_bend/r
    
#l=4 #length of s/c config gave 3.349
#d=2.5 #diameter
 #radius config gave 0.71
