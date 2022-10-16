#%%
# Import Dependencies
from IPython.display import display_markdown
from pyvlm import LatticeResult
from pyvlm import latticesystem_from_json

#%%
# Create Lattice System
jsonfilepath = '../files/Prandtl-D2.json'
lsys = latticesystem_from_json(jsonfilepath)
display_markdown(lsys)

#%%
# Design Point Result
alpha = 0.0
speed = 13.0
rho = 1.145

lres = LatticeResult('Design Point', lsys)
lres.set_density(rho=rho)
lres.set_state(alpha=alpha, speed=speed)

display_markdown(lres)
display_markdown(lres.surface_loads)
display_markdown(lres.stability_derivatives)

#%%
# Roll Case Result
alpha = 0.0
speed = 13.0
pbo2V = 0.01
rho = 1.145

lres = LatticeResult('Roll Case', lsys)
lres.set_density(rho=rho)
lres.set_state(alpha=alpha, speed=speed, pbo2V=pbo2V)

display_markdown(lres)
display_markdown(lres.surface_loads)
display_markdown(lres.stability_derivatives)

#%%
# Plots
axp = None
axp = lres.plot_phi_distribution(ax=axp)
_ = axp.set_ylabel('Phi Distribution')
_ = axp.set_xlabel('Span Position')

axl = None
axl = lres.plot_trefftz_lift_force_distribution(ax=axl)
axl = lres.plot_strip_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

axd = None
axd = lres.plot_trefftz_drag_force_distribution(ax=axd)
axd = lres.plot_strip_drag_force_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

axw = None
axw = lres.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%%
# Display Strip Geometry
display_markdown(lsys.strip_geometry)
