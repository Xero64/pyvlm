#%%
# Import Dependencies
from IPython.display import display_markdown

from pyvlm import LatticeOptimum, latticesystem_from_json
from pyvlm.tools import Bell

#%%
# Low AR Wing
jsonfilepath = '../files/Straight_Wing_Cosine_100.json'
lsys = latticesystem_from_json(jsonfilepath)
display_markdown(lsys)

#%%
# Bell
y = [lsys.srfcs[0].pnts[i][0].y for i in range(len(lsys.srfcs[0].pnts))]

bll = Bell(lsys.bref, y)
bll.set_lift(1.0)
bll.set_ym(lsys.srfcs[0].strpy)
l = bll.return_phi()

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
axl = bll.plot_lift_force_distribution(ax=axl)

axd = None
axd = lopt.plot_trefftz_drag_force_distribution(ax=axd)
axd = bll.plot_drag_force_distribution(ax=axd)

axw = None
axw = lopt.plot_trefftz_wash_distribution(ax=axw)
axw = bll.plot_trefftz_wash_distribution(ax=axw)
