#%% Import Dependencies
from IPython.display import display_markdown
from pyvlm_files import load_package_file
from pyvlm.outputs.msh import latticeresult_to_msh
from pyvlm.outputs.prf import latticeresult_to_prf

#%% Import Geometry
jsonfilepath = 'Test_aircraft.json'
lsys = load_package_file(jsonfilepath)

#%% Display System
display_markdown(lsys)

#%% Display Results

for case in lsys.results:
    lres = lsys.results[case]
    display_markdown(lres)

#%% Mesh File Output
lres = lsys.results['Positive 1g Cruise + 15deg Side Slip']
latticeresult_to_msh(lres, 'Test_aicraft.msh')

#%% Pessure File Output
latticeresult_to_prf(lsys, 'Test_aicraft_pressures.json')

#%% 5g Trim Case
ltrm = lsys.results['Positive 5g Dive']

#%% Plot Lift Distribution
axl = ltrm.plot_trefftz_lift_distribution()

#%% Plot Y Force Distribution
axy = ltrm.plot_trefftz_yforce_distribution()

#%% Print Strip Forces
display_markdown(ltrm.strip_forces)

#%% Print Strip Coefficients
display_markdown(ltrm.strip_coefficients)

#%% Print Panel Forces
display_markdown(ltrm.panel_forces)

#%% Print Total Loads
display_markdown(ltrm.surface_loads)
