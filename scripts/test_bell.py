#%% Import Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm_files import load_package_file
from pyvlm.tools import Bell

#%% Low AR Wing
jsonfilename = "Straight_Wing_Cosine_100.json"
lsys = load_package_file(jsonfilename)
print(lsys)

#%% Bell

# y = [strpi.pnti.y for strpi in lsys.strps]

y = [pnt.y for pnt in lsys.srfcs[0].pnts[:, 0].transpose().tolist()[0]]

bll = Bell(lsys.bref, y)
bll.set_lift(1.0)
bll.set_ym(lsys.srfcs[0].strpy)
l = bll.return_phi()

#%% Low AR Wing Optimum

lres = LatticeResult('Low AR Wing', lsys)
lres.set_state()
lres.set_lift_distribution(l, rho=1.0, speed=1.0)
print(lres)

#%% Plots

axl = None
axl = lres.plot_trefftz_lift_distribution(ax=axl)
axl = bll.plot_lift_distribution(ax=axl)

axd = None
axd = lres.plot_trefftz_drag_distribution(ax=axd)
axd = bll.plot_drag_distribution(ax=axd)

axw = None
axw = lres.plot_trefftz_wash_distribution(ax=axw)
axw = bll.plot_trefftz_wash_distribution(ax=axw)
