#%% Import Dependencies
from pyvlm import latticesystem_from_json, LatticeResult, LatticeOptimum
from pyvlm.tools import Elliptical

#%% Low AR Wing

jsonfilepath1 = r"packages\pyvlm\files\Straight_Wing_Equal_20.json"
lsys = latticesystem_from_json(jsonfilepath1)

print(f'bref = {lsys.bref:g}')
print(f'cref = {lsys.cref:g}')
print(f'sref = {lsys.sref:g}')
print(f'rref = {lsys.rref:.3f}')

#%% Elliptical

# y = [strpi.pnti.y for strpi in lsys.strps]

y = [pnt.y for pnt in lsys.srfcs[0].pnts[:, 0].transpose().tolist()[0]]

ell = Elliptical(lsys.bref, y)
ell.set_lift(1.0)

phi = ell.return_phi()

#%% Low AR Wing Optimum

lres = LatticeResult('Low AR Wing', lsys)
lres.set_conditions()
lres.set_phi(phi)
lres.print_aerodynamic_coefficients()

#%% Elliptical

ell = Elliptical(lsys.bref, lres.strpy)
ell.set_lift(1.0)

print(ell)

#%% Plots

axl = None
axl = lres.plot_trefftz_lift_distribution(ax=axl)
axl = ell.plot_lift_distribution(ax=axl)

axd = None
axd = lres.plot_trefftz_drag_distribution(ax=axd)
axd = ell.plot_drag_distribution(ax=axd)

axw = None
axw = lres.plot_trefftz_wash_distribution(ax=axw)
axw = ell.plot_wash_distribution(ax=axw)

#%% Variables

bvg = lsys.bvg.tolist()

#%%
