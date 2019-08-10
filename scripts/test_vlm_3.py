#%% Import Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm.files import load_package_file

#%% Low AR Wing

jsonfilename = "Test.json"
lsys = load_package_file(jsonfilename)
lsys_new = load_package_file(jsonfilename)

print(f'bref = {lsys.bref:g}')
print(f'cref = {lsys.cref:g}')
print(f'sref = {lsys.sref:g}')
print(f'rref = {lsys.rref:.3f}')

lres_org = LatticeResult('Initial', lsys)
lres_org.set_conditions(alpha=3.0)

Lspec = lres_org.CL*lres_org.qfs*lsys.sref

#%% Low AR Optimal Lift Distribution

lopt = LatticeOptimum('Optimal', lsys_new)
lopt.set_conditions()
lopt.add_constraint('L', Lspec)

phi, lam = lopt.optimum_lift_distribution()

# print(phi)

# phi, lam, Di, L, l = lsys.optimum_lift_distribution(Lspec)

# print(phi)

# print(f'phi1 = {phi1}')
# print(f'lam1 = {lam1}')
# print(f'Di = {Di}')
# print(f'L = {L}')
# print(f'l = {l}')

lres_opt = LatticeResult('Optimal', lsys)
lres_opt.set_conditions(alpha=3.0)
lres_opt.set_phi(phi)
lres_opt.print_aerodynamic_coefficients()

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
axa.plot(lres_opt.strpy, al1, label='al1')
axa.plot(lres_opt.strpy, al2, label='al2')
axa.legend()
