#%%
# Import Dependencies
from IPython.display import display_markdown
from pyvlm import LatticeOptimum, LatticeSystem
from pyvlm.tools import Elliptical

#%%
# Low AR Wing
jsonfilepath = '../files/Straight_Wing_Cosine_100.json'
lsys = LatticeSystem.from_json(jsonfilepath)
display_markdown(lsys)

#%%
# Elliptical
y = [lsys.srfcs[0].pnts[i][0].y for i in range(len(lsys.srfcs[0].pnts))]

ell = Elliptical(lsys.bref, y)
ell.set_lift(1.0)
ell.set_ym(lsys.srfcs[0].strpy)
l = ell.return_phi()

#%%
# Low AR Wing Optimum
lopt = LatticeOptimum('Low AR Wing', lsys)
lopt.set_state()
lopt.set_target_lift_force_distribution(l, rho=1.0, speed=1.0)
display_markdown(lopt)

#%%
# Plots
axl = None
axl = lopt.plot_trefftz_lift_force_distribution(ax=axl)
axl = ell.plot_lift_force_distribution(ax=axl)

axd = None
axd = lopt.plot_trefftz_drag_force_distribution(ax=axd)
axd = ell.plot_drag_force_distribution(ax=axd)

axw = None
axw = lopt.plot_trefftz_wash_distribution(ax=axw)
axw = ell.plot_trefftz_wash_distribution(ax=axw)
