#%% Load Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm.files import load_package_file

#%% Create Lattice System
jsonfilename = "Straight_Wing_Cosine_100.json"
lsys = load_package_file(jsonfilename)
print(lsys)

#%% Original Strip Geometry
print(lsys.strip_geometry)

#%% Original Strip Geometry
print(lsys.panel_geometry)

#%% Original Case

lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=3.0)#, pbo2V=0.002)
print(lres_org)

# lres_org.print_total_loads()
# lres_org.print_aerodynamic_coefficients()

#%% Original Strip Forces
print(lres_org.strip_forces)

#%% Original Strip Coefficients
print(lres_org.strip_coefficients)

#%% Original Panel Forces
print(lres_org.panel_forces)

#%% Plot Distribution
axl = lres_org.plot_trefftz_lift_distribution()

#%% Plot Wash
axw = lres_org.plot_trefftz_wash_distribution()

#%% Plot Near Field Velocity
axw = lres_org.plot_panel_near_field_velocities(component='z')

#%% Print Results
print(lres_org)

#%% Print Stability Derivatives
print(lres_org.stability_derivatives)

#%% Print Control Derivatives
print(lres_org.control_derivatives)
