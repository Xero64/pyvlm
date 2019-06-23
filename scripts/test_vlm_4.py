#%% Import Dependencies

from pyvlm import latticesystem_from_json, LatticeResult, LatticeOptimum
from pyvlm.tools import bell_lift_distribution

#%% Parameters
V = 18.0 # m/s
rho = 1.206 # kg/m**3
q = rho*V**2/2 # Pa
print(f'q = {q:.2f} Pa')

m = 6.577089 # kg
g = 9.80655 # m/s**2
W = m*g # N
print(f'W = {W:.2f} N')

#%% Low AR Wing

jsonfilepath = r"packages\pyvlm\files\Test_rhofw.json"
lsys = latticesystem_from_json(jsonfilepath)
lsys_opt = latticesystem_from_json(jsonfilepath)
lsys_bll = latticesystem_from_json(jsonfilepath)

print(f'bref = {lsys.bref:g}')
print(f'cref = {lsys.cref:g}')
print(f'sref = {lsys.sref:g}')
print(f'rref = {lsys.rref:.3f}')

lres_org = LatticeResult('Initial', lsys)
lres_org.set_conditions(speed=V, rho=rho)

l = bell_lift_distribution(lres_org.strpy, lsys.bref, W)

#%% Bell Shaped Lift Distribution

CL = W/q/lsys.sref

print(f'CL = {CL:.5f}')

lopt_bll = LatticeOptimum('Bell', lsys_bll)
lopt_bll.set_conditions(speed=V, rho=rho)
lopt_bll.set_lift_distribution(l, rho, V)
lopt_bll.add_record('l', strplst='Mirrored')
lopt_bll.print_report()

lres_bll = LatticeResult('Bell', lsys)
lres_bll.set_conditions(speed=V, rho=rho)
lres_bll.set_lift_distribution(l, rho, V)
lres_bll.print_aerodynamic_coefficients()

#%% Optimal Lift Distribution

lopt = LatticeOptimum('Optimal', lsys_opt)
lopt.set_conditions(speed=V, rho=rho)
lopt.add_constraint('L', W)
lopt.add_constraint('l', 20.566097, strplst='Mirrored')
# lopt.add_record('l', strplst='Mirrored')
lopt.optimum_lift_distribution()
lopt.print_report()

lres_opt = LatticeResult('Optimal', lsys)
lres_opt.set_conditions(speed=V, rho=rho)
lres_opt.set_phi(lopt.phi)
lres_opt.print_aerodynamic_coefficients()

#%% Plots

axp = None
axp = lres_org.plot_phi_distribution(ax=axp)
axp = lres_bll.plot_phi_distribution(ax=axp)
axp = lres_opt.plot_phi_distribution(ax=axp)

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_bll.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt.plot_trefftz_lift_distribution(ax=axl)

axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lres_bll.plot_trefftz_drag_distribution(ax=axd)
axd = lres_opt.plot_trefftz_drag_distribution(ax=axd)

axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lres_bll.plot_trefftz_wash_distribution(ax=axw)
axw = lres_opt.plot_trefftz_wash_distribution(ax=axw)

#%% Optimal Strip Twist

al_opt = lopt.optimum_strip_twist(crit=1e-1)
al_bll = lopt_bll.optimum_strip_twist(crit=1e-1)

#%% Specified String Twist

alspec = [8.3274, 8.5524, 8.7259, 8.8441, 8.9030, 8.8984, 8.8257, 8.6801, 8.4565, 8.1492, 7.7522,
          7.2592, 6.6634, 5.9579, 5.1362, 4.1927, 3.1253, 1.9394, 0.6589, -0.6417, -1.6726]

numal = len(alspec)

yspec = [i*lsys.bref/2/(numal-1) for i in range(numal)]

#%% Plots

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_bll.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt.plot_trefftz_lift_distribution(ax=axl)
axl = lopt.res.plot_trefftz_lift_distribution(ax=axl)

axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lres_bll.plot_trefftz_drag_distribution(ax=axd)
axd = lres_opt.plot_trefftz_drag_distribution(ax=axd)
axd = lopt.res.plot_trefftz_drag_distribution(ax=axd)

axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lres_bll.plot_trefftz_wash_distribution(ax=axw)
axw = lres_opt.plot_trefftz_wash_distribution(ax=axw)
axw = lopt.res.plot_trefftz_wash_distribution(ax=axw)

axa = lopt.plot_strip_twist_distribution()
axa.plot(lres_bll.strpy, al_bll, label='alpha Bell')
axa.plot(lres_opt.strpy, al_opt, label='alpha Optimum')
axa.plot(yspec, alspec, label='alpha Specified')
leg = axa.legend()
