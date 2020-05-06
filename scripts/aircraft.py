#%% Import Dependencies
from IPython.display import display_markdown
from py2md.classes import MDMatrix
from pyvlm import latticesystem_from_json
from pyvlm.outputs.msh import latticeresult_to_msh
from pyvlm.outputs.prf import latticeresult_to_prf
from pyvlm.classes import LatticeTrim
from numpy import zeros
from numpy.linalg import eig, solve
from math import cos, radians
from pyvlm.tools import StabilityApproximation

#%% Import Geometry
jsonfilepath = r'..\files\Aircraft.json'
lsys = latticesystem_from_json(jsonfilepath)

#%% Display System
display_markdown(lsys)

#%% Display Results
for case in lsys.results:
    lres = lsys.results[case]
    display_markdown(lres)

#%% Mesh File Output
lres = lsys.results['Positive 1g Cruise + 15deg Side Slip']
latticeresult_to_msh(lres, r'..\results\Aircraft.msh')

#%% Pessure File Output
latticeresult_to_prf(lsys, r'..\results\Aircraft_pressures.json')

#%% 5g Trim Case
ltrm = lsys.results['Positive 5g Dive']

#%% Plot Lift Distribution
axl = ltrm.plot_trefftz_lift_distribution()
axl = ltrm.plot_strip_lift_distribution(ax=axl)

#%% Plot Y Force Distribution
axy = ltrm.plot_trefftz_yforce_distribution()
axy = ltrm.plot_strip_yforce_distribution(ax=axy)

#%% Plot Drag Distribution
axd = ltrm.plot_trefftz_drag_distribution()
axd = ltrm.plot_strip_drag_distribution(ax=axd)

#%% Print Strip Forces
display_markdown(ltrm.strip_forces)

#%% Print Strip Coefficients
display_markdown(ltrm.strip_coefficients)

#%% Print Panel Forces
display_markdown(ltrm.panel_forces)

#%% Print Total Loads
display_markdown(ltrm.surface_loads)

#%% Trim CL to 0.8
CLt = 0.8
CYt = 0.0

ltrm2 = LatticeTrim(f'CL = {CLt}, CY = {CYt}', lsys)
ltrm2.set_targets(CLt = CLt, CYt = CYt)
ltrm2.set_trim_loads(trmmom=False)
ltrm2.trim()

display_markdown(ltrm2)

#%% Neutral Point
for case in lsys.results:
    lres = lsys.results[case]
    print(f'Case: {case:s}')
    print(f'Neutral Point: {lres.stres.neutral_point():.3f}')

display_markdown(lres.stability_derivatives)

#%% Centre of Pressure
for case in lsys.results:
    lres = lsys.results[case]
    print(f'Case: {case:s}')
    x, y, z = lres.nfres.centre_of_pressure('zx')
    print(f'Centre of Pressure: {x:.3f}, {y:.3f}, {z:.3f}')

#%% Neutral Point 2
for case in lsys.results:
    lres = lsys.results[case]
    print(f'Case: {case:s}')
    x, y, z = lres.stres.alpha.centre_of_pressure('zx')
    print(f'Neutral Point: {x:.3f}, {y:.3f}, {z:.3f}')

#%% State Derivative
# lsys.masses['Nominal CG'].Ixx = 500.0
# lsys.masses['Nominal CG'].Iyy = 500.0
# lsys.masses['Nominal CG'].Izz = 500.0

# A = lres.stres.system_aerodynamic_matrix()
# mat = MDMatrix('A', A, '.3f')
# display_markdown(mat)

# M = zeros((6, 6))
# M[0:3, 0:3] = lsys.masses['Nominal CG'].mass_matrix()
# M[3:6, 3:6] = lsys.masses['Nominal CG'].moment_of_inertia_matrix()

# mat = MDMatrix('M', M, '.3f')
# display_markdown(mat)

# newA = solve(M, A)
# mat = MDMatrix('newA', newA, '.3f')
# display_markdown(mat)

# lonA = zeros((4, 4))
# lonA[:3, :3] = newA[0::2, 0::2]
# lonA[0, 2] = 0.0
# lonA[1, 2] = lres.vfs.x
# lonA[0, 3] = -9.8065
# lonA[3, 2] = 1.0

# mat = MDMatrix('lonA', lonA, '.3f')
# display_markdown(mat)

# wlonA, vlonA = eig(lonA)

# wlonA = wlonA.reshape((4, 1))

# mat = MDMatrix('wlonA', wlonA, '.6f')
# display_markdown(mat)

# mat = MDMatrix('vlonA', vlonA, '.3f')
# display_markdown(mat)

# latA = zeros((4, 4))
# latA[:3, :3] = newA[1::2, 1::2]
# latA[0, 2] = latA[0, 2]-lres.vfs.x
# latA[0, 3] = -9.8065*cos(radians(lres.alpha))
# latA[3, 1] = 1.0

# mat = MDMatrix('latA', latA, '.3f')
# display_markdown(mat)

# wlatA, vlatA = eig(latA)

# wlatA = wlatA.reshape((4, 1))

# mat = MDMatrix('wlatA', wlatA, '.6f')
# display_markdown(mat)

# mat = MDMatrix('vlatA', vlatA, '.3f')
# display_markdown(mat)

# #%% Stability Approximations

# stabapprox = StabilityApproximation(lres, lsys.masses['Nominal CG'])

# print(stabapprox.phugoid())
# print(stabapprox.short_period())
# print(stabapprox.roll_subsidence())
# print(stabapprox.spiral())
# print(stabapprox.dutch_roll())
