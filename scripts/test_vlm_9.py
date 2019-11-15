#%% Load Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm.files import load_package_file

#%% Create Lattice System
jsonfilename = "Straight_Wing_Control.json"
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

#%% Plot Aileron Drag Derivative
ctres = lres_org.ctresp['aileron']

#%%
from pygeom.geom3d import Vector
from matplotlib.pyplot import figure

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
            force += ctres.nffrc[i, 0]
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

lres_ail = LatticeResult('Baseline', lsys)
lres_ail.set_state(alpha=3.0)
lres_ail.set_controls(aileron=2.0)
print(lres_ail)

axl = lres_ail.plot_trefftz_lift_distribution()
