#%% Load Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm.files import load_package_file

#%% Create Lattice System
jsonfilename = "Straight_Wing_Cosine_100.json"
lsys = load_package_file(jsonfilename)
print(lsys)

#%% Original Strip Geometry
lsys.print_strip_geometry()

#%% Original Strip Geometry
lsys.print_panel_geometry()

#%% Original Case

alpha = 3.0 # degrees
beta = 0.0 # degrees

lres_org = LatticeResult('Baseline', lsys)
lres_org.set_conditions(alpha=alpha, beta=beta, omx=0.002)
print(lres_org)

# lres_org.print_total_loads()
# lres_org.print_aerodynamic_coefficients()

#%% Original Strip Forces
lres_org.print_strip_forces()

#%% Original Strip Coefficients
lres_org.print_strip_coefficients()

#%% Original Panel Forces
lres_org.print_panel_forces()

#%% Plot Distribution

axl = lres_org.plot_trefftz_lift_distribution()

#%% Plot Wash

axw = lres_org.plot_trefftz_wash_distribution()

#%%

axw = lres_org.plot_panel_near_field_velocities(component='z')

#%%
