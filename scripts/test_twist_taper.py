#%% Load Dependencies
from IPython.display import display_markdown
from pyvlm import LatticeResult
from pyvlm import latticesystem_from_json

#%% Create Lattice System
jsonfilepath = '../files/Test_twist_taper.json'
lsys = latticesystem_from_json(jsonfilepath)
display_markdown(lsys)

#%% Original Strip Geometry
display_markdown(lsys.strip_geometry)

#%% Original Strip Geometry
display_markdown(lsys.panel_geometry)

#%% Original Case
lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=3.0)

#%% Original Strip Forces
display_markdown(lres_org.strip_forces)

#%% Original Strip Coefficients
display_markdown(lres_org.strip_coefficients)

#%% Original Panel Forces
display_markdown(lres_org.panel_forces)

#%% Print Result
display_markdown(lres_org)

#%% Plot Lift Distribution
axl = None
axl = lres_org.plot_trefftz_lift_force_distribution(ax=axl)
axl = lres_org.plot_strip_lift_force_distribution(ax=axl)

#%% Plot Drag Distribution
axd = None
axd = lres_org.plot_trefftz_drag_force_distribution(ax=axd)
axd = lres_org.plot_strip_drag_force_distribution(ax=axd)

#%% Plot Wash Distribution
axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
