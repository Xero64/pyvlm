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

lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=3.0)#, pbo2V=0.002)
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

#%% Plot Near Field Velocity
axw = lres_org.plot_panel_near_field_velocities(component='z')

#%% Print Results
print(lres_org)

#%% Print Stability Derivatives
lres_org.print_stability_derivatives()

#%% Print Control Derivatives
lres_org.print_control_derivatives()

#%% Checks
# print(lres_org.ofs.x*lsys.bref/2/lres_org.speed)
# print(lres_org.ofs.y*lsys.cref/2/lres_org.speed)
# print(lres_org.ofs.z*lsys.bref/2/lres_org.speed)