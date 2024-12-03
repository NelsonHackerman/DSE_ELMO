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

m=6519.77 #(new estimate)
l=7 #length of s/c
d=2 #diameter
r=d/2 #radius
l_ax_com=6*g*SF*m #axial load for compression
l_ax_ten=2.5*g*SF*m #axial load for tension
l_lat=2*g*SF*m #lateral load
m_bend=l/2*l_lat #axial bending moment from lateral loads
l_eq_ten=l_ax_ten+2*m_bend/r #axial tension load with bending
l_eq_com=l_ax_com+2*m_bend/r #axial compressive load with bending
f_ax=25 #axial natural frequency requirement
f_lat=10 #longitudinal natural frequency requirement