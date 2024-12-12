#material: AL 7075
E=71*10**9#N/m2
v=0.33
rho=2.8*10**3#kg/m3
Sig_tu=524*10**6#N/m2 ultimate tensile strength
Sig_ty=448*10**6#N/m2 yield tensile strength
g=9.81
SF=1.1 #table 11-54 SMAD
#falcon 9: axial -2 to 6g, long -2 to 2g + for compress
#arianes: axial -6 to 2.5g, long -1.8 to 1.8g - for compress
#overall axial -6 to 2.5g, long -2 to 2g (I'M USING + FOR TENSION)
#falcon 9: axial >25Hz, lat >10Hz
#arianes: axial >20Hz, lat >6Hz
#overall axial >25Hz, lat >10Hz

option='d'
if option=='a':
    morb=2778.26
    mtot=8629.51
if option=='b':
    m=13012.01
if option=='c':
    morb=4618.22
    mtot=13618.38
if option=='d':
    m=18762.26
    
l=4 #length of s/c config gave 3.349
d=2.5 #diameter
r=d/2 #radius config gave 0.71
l_ax_com=6*g*SF #axial load for compression
l_ax_ten=2.5*g*SF #axial load for tension
l_lat=2*g*SF #lateral load
m_bend=l/2*l_lat #axial bending moment from lateral loads
l_eq_ten=l_ax_ten+2*m_bend/r #axial tension load with bending
l_eq_com=l_ax_com+2*m_bend/r #axial compressive load with bending
f_ax=25 #axial natural frequency requirement
f_lat=10 #longitudinal natural frequency requirement