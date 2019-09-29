#%% Import Dependencies
from pyvlm.tools import full_cosine_spacing
from pyvlm.tools.camber import NACA4, NACA6Series

#%% Calculate Spacing
naca6 = NACA6Series('65(1)212')

xc = full_cosine_spacing(8*4+2)
sl = [naca6.return_camber_slope(xci) for xci in xc]

for i in range(len(xc)):
    xci = xc[i]
    sli = sl[i]
    print(f'{xci:.6f}, {sli:.6f}')
