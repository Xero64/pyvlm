#%% Import Dependencies
from IPython.display import display_markdown
from pyvlm import latticesystem_from_json
from pyvlm.outputs.msh import latticeresult_to_msh
from pyvlm.outputs.prf import latticeresult_to_prf
from pyvlm.classes import LatticeTrim

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
