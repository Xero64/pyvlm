#%% Load Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm.files import load_package_file

#%% Create Lattice System
jsonfilename = "Test_twist.json"
lsys = load_package_file(jsonfilename)
print(lsys)

#%% Original Strip Geometry
lsys.print_strip_geometry()

#%% Original Strip Geometry
lsys.print_panel_geometry()

#%% Original Case
lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=3.0)

#%% Print Result
print(lres_org)

#%% Original Strip Forces
lres_org.print_strip_forces()

#%% Original Strip Coefficients
lres_org.print_strip_coefficients()

#%% Original Panel Forces
lres_org.print_panel_forces()
