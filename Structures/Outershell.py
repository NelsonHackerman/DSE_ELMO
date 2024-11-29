import numpy as np


#option 1: monocoque- no ring or longitudinal stiffeners
#material: AL 7075
E=71*10**9#N/m2
v=0.33
rho=2.8*10**3#kg/m3
F_tu=524*10**6#N/m2 ultimate tensile strength
F_ty=448*10**6#N/m2 yield tensile strength
#falcon 9: axial -2 to 6g, long -2 to 2g
#arianes: axial -6 to 2.5g, long -1.8 to 1.8g
#overall axial -6 to 6g, long -2 to 2g
#falcon 9: axial >25Hz, lat >10Hz
#arianes: axial >20Hz, lat >6Hz
#overall axial >25Hz, lat >10Hz
f_ax=25
f_lat=10
l=np.sqrt(E/rho)*0.25/f_ax
print('Required length for natural frequency monocoque : ',l)
