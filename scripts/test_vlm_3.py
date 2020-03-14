#%% Import Dependencies
from IPython.display import display
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm import latticesystem_from_json

#%% Low AR Wing
jsonfilepath = r'..\files\Test.json'
lsys = latticesystem_from_json(jsonfilepath)
display(lsys)

lsys_opt = latticesystem_from_json(jsonfilepath)

lres_org = LatticeResult('Initial', lsys)
lres_org.set_state(alpha=3.0)
display(lres_org)

Lspec = lres_org.trres.CL*lres_org.qfs*lsys.sref

#%% Low AR Optimal Lift Distribution
lopt = LatticeOptimum('Optimal', lsys_opt)
lopt.set_state()
lopt.add_constraint('L', Lspec)
phi, lam = lopt.optimum_lift_distribution()
display(lopt)

#%% Original Lift Distribution Plot
axl = lres_org.plot_trefftz_lift_distribution()
axl = lopt.plot_trefftz_lift_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')
_ = axl.legend()

#%% Original Drag Distribution Plot
axd = lres_org.plot_trefftz_drag_distribution()
axd = lopt.plot_trefftz_drag_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')
_ = axd.legend()

#%% Original Wash Distribution Plot
axw = lres_org.plot_trefftz_wash_distribution()
axw = lopt.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')
_ = axw.legend()

#%% Optimal Strip Twist
al1 = lopt.optimum_strip_twist(crit=1e-1)
al2 = lopt.optimum_strip_twist(crit=1e-2)

#%% Lift Distribution Plot
axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lopt.plot_trefftz_lift_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%% Drag Distribution Plot
axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lopt.plot_trefftz_drag_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

#%% Wash Distribution Plot
axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lopt.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%% Geometric Twist Distribution Plot
axa = lopt.plot_strip_twist_distribution()
axa.plot(lsys.srfcs[0].strpy, al1, label='al1')
axa.plot(lsys.srfcs[0].strpy, al2, label='al2')
_ = axw.set_ylabel('Strip Twist Distribution [deg]')
_ = axw.set_xlabel('Span Position')
_ = axa.legend()

#%% New 0deg Results
lres_0deg = LatticeResult('0deg Result', lsys_opt)
lres_0deg.set_state(alpha=0.0, speed=1.0)
display(lres_0deg)

#%% New 6deg Results
lres_6deg = LatticeResult('6deg Result', lsys_opt)
lres_6deg.set_state(alpha=6.0, speed=1.0)
display(lres_6deg)

#%% New Lift Distribution Plot
axl = None
axl = lres_0deg.plot_trefftz_lift_distribution(ax=axl)
axl = lres_6deg.plot_trefftz_lift_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%% New Drag Distribution Plot
axd = None
axd = lres_0deg.plot_trefftz_drag_distribution(ax=axd)
axd = lres_6deg.plot_trefftz_drag_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

#%% New Wash Distribution Plot
axw = None
axw = lres_0deg.plot_trefftz_wash_distribution(ax=axw)
axw = lres_6deg.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')
