#%%
# Load Dependencies
from pyvlm import LatticeResult, latticesystem_from_json

#%%
# Create Lattice System
jsonfilepath = '../files/Test_twist_sweep_taper.json'
lsys = latticesystem_from_json(jsonfilepath)
print(lsys)

#%%
# Original Strip Geometry
print(lsys.strip_geometry)

#%%
# Original Strip Geometry
print(lsys.panel_geometry)

#%%
# Original Case
lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=3.0)

#%%
# Original Strip Forces
print(lres_org.strip_forces)

#%%
# Original Strip Coefficients
print(lres_org.strip_coefficients)

#%%
# Original Panel Forces
print(lres_org.panel_forces)

#%%
# Print Result
print(lres_org)
