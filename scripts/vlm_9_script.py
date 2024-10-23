#%%
# Load Dependencies
from IPython.display import display

from pyvlm import LatticeResult, latticesystem_from_json

#%%
# Create Lattice System
jsonfilepath = '../files/Straight_Wing_Control.json'
lsys = latticesystem_from_json(jsonfilepath)
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
lres_org.set_state(alpha=3.0)#, pbo2V=0.002)
display(lres_org)

# lres_org.print_total_loads()
# lres_org.print_aerodynamic_coefficients()

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
# Print Results
display(lres_org)

#%%
# Print Stability Derivatives
display(lres_org.stability_derivatives)

#%%
# Print Control Derivatives
display(lres_org.control_derivatives)

#%%
# Plot Aileron Drag Derivative
ctres = lres_org.ctresp['aileron']

from matplotlib.pyplot import figure
#%%
# Plot Aileron Deflection Plots
from pygeom.geom3d import Vector

figl = figure(figsize=(12, 8))
axl = figl.gca()
axl.grid(True)

figd = figure(figsize=(12, 8))
axd = figd.gca()
axd.grid(True)

for srfc in ctres.res.sys.srfcs:
    y = []
    d = []
    l = []
    for strp in srfc.strps:
        y.append(strp.pnti.y)
        force = Vector(0.0, 0.0, 0.0)
        for pnl in strp.pnls:
            i = pnl.lpid
            force += ctres.nffrc[i]
        locfrc = ctres.res.acs.vector_to_local(force)
        li = locfrc.z
        l.append(li)
        di = locfrc.x
        d.append(di)
    if len(l) > 0:
        label = ctres.res.name + ' for ' + srfc.name
        axl.plot(y, l, label=label)
    if len(d) > 0:
        label = ctres.res.name + ' for ' + srfc.name
        axd.plot(y, d, label=label)

_ = axl.legend()
_ = axd.legend()

#%%
# Aileron Deflection

lres_ail = LatticeResult('Baseline', lsys)
lres_ail.set_state(alpha=3.0)
lres_ail.set_controls(aileron=2.0)
display(lres_ail)

axl = lres_ail.plot_trefftz_lift_force_distribution()
