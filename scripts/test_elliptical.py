#%% Import Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm.files import load_package_file
from pyvlm.tools import Elliptical

#%% Low AR Wing
jsonfilename = "Straight_Wing_Cosine_100.json"
lsys = load_package_file(jsonfilename)
print(lsys)

#%% Elliptical

# y = [strpi.pnti.y for strpi in lsys.strps]

y = [pnt.y for pnt in lsys.srfcs[0].pnts[:, 0].transpose().tolist()[0]]

ell = Elliptical(lsys.bref, y)
ell.set_lift(1.0)
ell.set_ym(lsys.strpy)
l = ell.return_phi()

#%% Low AR Wing Optimum

lres = LatticeResult('Low AR Wing', lsys)
lres.set_conditions()
lres.set_lift_distribution(l, rho=1.0, speed=1.0)
print(lres)

#%% Plots

axl = None
axl = lres.plot_trefftz_lift_distribution(ax=axl)
axl = ell.plot_lift_distribution(ax=axl)

axd = None
axd = lres.plot_trefftz_drag_distribution(ax=axd)
axd = ell.plot_drag_distribution(ax=axd)

axw1 = None
axw1 = lres.plot_trefftz_wash_distribution(ax=axw1)
axw2 = None
axw2 = ell.plot_trefftz_wash_distribution(ax=axw2)

#%% Variables

well = ell.trefftz_wash_distribution()
wres = lres.trwsh

from matplotlib.pyplot import figure

fig = figure()
ax = fig.gca()
ax.grid(True)
ax.plot([0.0], [0.0])
ax.plot(ell.y, well, label='Theory')
ax.plot(lsys.strpy, wres, label='Result')
l = ax.legend()

# # bvg = lsys.bvg.tolist()


#%%
