#%% Import Dependencies
from pyvlm import LatticeOptimum
from pyvlm import latticesystem_from_json
from pyvlm.tools import Elliptical

#%% Low AR Wing
jsonfilepath = '../files/Straight_Wing_Cosine_100.json'
lsys = latticesystem_from_json(jsonfilepath)
print(lsys)

#%% Elliptical

# y = [strpi.pnti.y for strpi in lsys.strps]

y = [pnt.y for pnt in lsys.srfcs[0].pnts[:, 0].transpose().tolist()[0]]

ell = Elliptical(lsys.bref, y)
ell.set_lift(1.0)
ell.set_ym(lsys.srfcs[0].strpy)
l = ell.return_phi()

#%% Low AR Wing Optimum

lopt = LatticeOptimum('Low AR Wing', lsys)
lopt.set_state()
lopt.set_target_lift_force_distribution(l, rho=1.0, speed=1.0)
print(lopt)

#%% Plots

axl = None
axl = lopt.plot_trefftz_lift_force_distribution(ax=axl)
axl = ell.plot_lift_force_distribution(ax=axl)

axd = None
axd = lopt.plot_trefftz_drag_force_distribution(ax=axd)
axd = ell.plot_drag_force_distribution(ax=axd)

axw1 = None
axw1 = lopt.plot_trefftz_wash_distribution(ax=axw1)
axw2 = None
axw2 = ell.plot_trefftz_wash_distribution(ax=axw2)

#%% Variables

well = ell.trefftz_wash_distribution()
wres = lopt.trres.trwsh

from matplotlib.pyplot import figure

fig = figure()
ax = fig.gca()
ax.grid(True)
ax.plot([0.0], [0.0])
ax.plot(ell.y, well, label='Theory')
ax.plot(lsys.srfcs[0].strpy, wres, label='Result')
l = ax.legend()

# # bvg = lsys.bvg.tolist()


#%%
