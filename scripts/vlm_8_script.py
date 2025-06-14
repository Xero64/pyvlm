#%%
# Load Dependencies
from IPython.display import display
from pyvlm import LatticeResult, LatticeSystem

#%%
# Create Lattice System
jsonfilepath = '../files/Straight_Wing_Cosine_100.json'
lsys = LatticeSystem.from_json(jsonfilepath)
display(lsys)

#%%
# Original Strip Geometry
display(lsys.strip_geometry)

#%%
# Original Strip Geometry
display(lsys.panel_geometry)

#%%
# Original Case
lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=3.0)#, pbo2v=0.002)
display(lres_org)

#%%
# Original Strip Forces
display(lres_org.strip_forces)

#%%
# Original Strip Coefficients
display(lres_org.strip_coefficients)

#%%
# Original Panel Forces
display(lres_org.panel_forces)

#%%
# Plot Distribution
axl = lres_org.plot_trefftz_lift_force_distribution()

#%%
# Plot Wash
axw = lres_org.plot_trefftz_wash_distribution()

#%%
# Plot Near Field Velocity
axw = lres_org.plot_panel_near_field_velocities(component='z')

#%%
# Print Results
display(lres_org)

#%%
# Print Stability Derivatives
display(lres_org.stability_derivatives)
