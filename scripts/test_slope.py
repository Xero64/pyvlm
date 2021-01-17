#%% Import Dependencies
from pyvlm.tools import full_cosine_spacing
from pyvlm.tools.camber import NACA4
from pyvlm.tools.airfoil import airfoil_from_dat

#%% Calculate Spacing
naca4 = NACA4('2412')

xc = full_cosine_spacing(8*4+2)
sl = [naca4.return_camber_slope(xci) for xci in xc]

for xci, sli in zip(xc, sl):
    print(f'{xci:.6f}, {sli:.6f}')

airfoil = airfoil_from_dat('../files/rhofw_root.dat')
sla = [airfoil.return_camber_slope(xci) for xci in xc]

print()
print(airfoil.name)
for xci, sli in zip(xc, sla):
    print(f'{xci:.6f}, {sli:.6f}')
