#%% Import Dependencies
from IPython.display import display
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm import latticesystem_from_json

#%% Low AR Wing
jsonfilepath = r'..\files\Test.json'
lsys = latticesystem_from_json(jsonfilepath)
display(lsys)

lsys_new = latticesystem_from_json(jsonfilepath)

lres_org = LatticeResult('Initial', lsys)
lres_org.set_state(alpha=3.0)
display(lres_org)

Lspec = lres_org.trres.CL*lres_org.qfs*lsys.sref

#%% Low AR Optimal Lift Distribution

lopt = LatticeOptimum('Optimal', lsys_new)
lopt.set_conditions()
lopt.add_constraint('L', Lspec)
phi, lam = lopt.optimum_lift_distribution()
display(lopt)

lres_opt = LatticeResult('Optimal', lsys)
lres_opt.set_state(alpha=3.0)
lres_opt.set_phi(phi)
display(lres_opt)

#%% Plots

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt.plot_trefftz_lift_distribution(ax=axl)

axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lres_opt.plot_trefftz_drag_distribution(ax=axd)

axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lres_opt.plot_trefftz_wash_distribution(ax=axw)

#%% Optimal Strip Twist

al1 = lopt.optimum_strip_twist(crit=1e-1)
al2 = lopt.optimum_strip_twist(crit=1e-2)

#%% Plots

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt.plot_trefftz_lift_distribution(ax=axl)
axl = lopt.res.plot_trefftz_lift_distribution(ax=axl)

axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lres_opt.plot_trefftz_drag_distribution(ax=axd)
axd = lopt.res.plot_trefftz_drag_distribution(ax=axd)

axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lres_opt.plot_trefftz_wash_distribution(ax=axw)
axw = lopt.res.plot_trefftz_wash_distribution(ax=axw)

axa = lopt.plot_strip_twist_distribution()
axa.plot(lsys.srfcs[0].strpy, al1, label='al1')
axa.plot(lsys.srfcs[0].strpy, al2, label='al2')
axa.legend()

#%% New Results

lres_0deg = LatticeResult('0deg Result', lsys_new)
lres_0deg.set_state(alpha=0.0, speed=1.0)
display(lres_0deg)

lres_6deg = LatticeResult('6deg Result', lsys_new)
lres_6deg.set_state(alpha=6.0, speed=1.0)
display(lres_6deg)

#%% New Plots

axl = None
axl = lres_0deg.plot_trefftz_lift_distribution(ax=axl)
axl = lres_6deg.plot_trefftz_lift_distribution(ax=axl)

axd = None
axd = lres_0deg.plot_trefftz_drag_distribution(ax=axd)
axd = lres_6deg.plot_trefftz_drag_distribution(ax=axd)

axw = None
axw = lres_0deg.plot_trefftz_wash_distribution(ax=axw)
axw = lres_6deg.plot_trefftz_wash_distribution(ax=axw)
