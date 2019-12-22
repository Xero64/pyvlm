#%% Import Dependencies
from IPython.display import display_markdown
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm_files import load_package_file
from pyvlm.tools import bell_lift_distribution

#%% Parameters
m = 6.577089 # kg
g = 9.80655 # m/s**2
W = m*g # N
print(f'W = {W:.2f} N')

S = 0.9375 # m**2
print(f'S = {S:.4f} m**2')

CL = 0.6
print(f'CL = {CL:.1f}')
rho = 1.206 # kg/m**3
print(f'rho = {rho:.3f} kg/m**3')
# rho = 1.0
V = (W/S/rho/CL*2)**0.5
print(f'V = {V:.4f} m/s')

q = rho*V**2/2 # Pa
print(f'q = {q:.2f} Pa')

#%% Low AR Wing

jsonfilename = "Test_rhofw.json"
lsys = load_package_file(jsonfilename)
lsys_opt = load_package_file(jsonfilename)
lsys_bll = load_package_file(jsonfilename)
display_markdown(lsys)

lres_org = LatticeResult('Initial', lsys)
lres_org.set_state(speed=V)
lres_org.set_density(rho=rho)
display_markdown(lres_org)

l = bell_lift_distribution(lsys.srfcs[0].strpy, lsys.bref, W)

#%% Bell Shaped Lift Distribution

lopt_bll = LatticeOptimum('Bell', lsys_bll)
lopt_bll.set_conditions(speed=V, rho=rho)
lopt_bll.set_lift_distribution(l, rho, V)
lopt_bll.add_record('l', strplst='Mirrored')
display_markdown(lopt_bll)

lres_bll = LatticeResult('Bell', lsys)
lres_bll.set_state(speed=V)
lres_bll.set_density(rho=rho)
lres_bll.set_lift_distribution(l, rho, V)
display_markdown(lres_bll)

#%% Optimal Lift Distribution

lopt = LatticeOptimum('Optimal', lsys_opt)
lopt.set_conditions(speed=V, rho=rho)
lopt.add_constraint('L', W)
lopt.add_constraint('l', 20.566097, strplst='Mirrored')
# lopt.add_record('l', strplst='Mirrored')
lopt.optimum_lift_distribution()
display_markdown(lopt)

lres_opt = LatticeResult('Optimal', lsys)
lres_opt.set_state(speed=V)
lres_opt.set_density(rho=rho)
lres_opt.set_phi(lopt.phi)
display_markdown(lres_opt)

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

almirr = []
ymirr = []

for i in range(numal-1, 0, -1):
    almirr.append(alspec[i])
    ymirr.append(-yspec[i])

alspec = almirr+alspec
yspec = ymirr+yspec

#%% Bell Downwash

w = [W/rho/V/lsys_bll.bref*3/2*((2*yi/lsys_bll.bref)**2-0.5) for yi in yspec]

print(min(w)/min(lres_bll.trres.trwsh))

#%% Plots

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_bll.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt.plot_trefftz_lift_distribution(ax=axl)
axl = lopt.res.plot_trefftz_lift_distribution(ax=axl)
axl = lopt.res.plot_strip_lift_distribution(ax=axl)

axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lres_bll.plot_trefftz_drag_distribution(ax=axd)
axd = lres_opt.plot_trefftz_drag_distribution(ax=axd)
axd = lopt.res.plot_trefftz_drag_distribution(ax=axd)
axd = lopt.res.plot_strip_drag_distribution(ax=axd)

axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lres_bll.plot_trefftz_wash_distribution(ax=axw)
axw = lres_opt.plot_trefftz_wash_distribution(ax=axw)
axw = lopt.res.plot_trefftz_wash_distribution(ax=axw)
axw.plot(yspec, w, label='Bell Wash')
axw.legend()

axa = lopt.plot_strip_twist_distribution()
axa.plot(lsys.srfcs[0].strpy, al_bll, label='alpha Bell')
axa.plot(lsys.srfcs[0].strpy, al_opt, label='alpha Optimum')
axa.plot(yspec, alspec, label='alpha Specified')
leg = axa.legend()

#%% Induced Angle

from math import atan2, degrees, radians, pi
from matplotlib.pyplot import figure
from pymath.function import Function

alf = Function(yspec, alspec)
als = alf.linear_interp(lsys.srfcs[0].strpy)

ali = [degrees(atan2(lres_bll.trres.trwsh[i], V)) for i in range(len(lres_bll.trres.trwsh))]
al0 = [1.0-abs(yi)*2/lsys_bll.bref for yi in lsys.srfcs[0].strpy]

alp = [als[i]+ali[i]+al0[i] for i in range(len(lsys.srfcs[0].strpy))]

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(lsys.srfcs[0].strpy, ali, label='Induced Angle')
ax.plot(lsys.srfcs[0].strpy, al0, label='Zero Lift Angle')
ax.plot(lsys.srfcs[0].strpy, als, label='Specified Angle')
ax.plot(lsys.srfcs[0].strpy, alp, label='Total Angle')
leg = ax.legend()

#%% Local Lift

cmax = 0.40005
cmin = 0.10008

c = [cmax-abs(yi)*2/lsys_bll.bref*(cmax-cmin) for yi in lsys.srfcs[0].strpy]

cla = 2*pi*0.91

l = [cla*c[i]*radians(alp[i])*q for i in range(len(lsys.srfcs[0].strpy))]

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(lsys.srfcs[0].strpy, l, label='Lift Distribution')
ax = lres_bll.plot_trefftz_lift_distribution(ax=ax)
leg = ax.legend()

#%% Plot Specified Twist

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(yspec, alspec)
ax.set_title('Wing Geometric Twist as a function of Span')
ax.set_xlabel('Span Coordinate - y - [m]')
_ = ax.set_ylabel('Geometric Twist Angle - $\\alpha_g$ - [deg]')

#%% Sum of Drag

stdrg = 0.0
for i in range(lopt.res.stripres.drag.shape[0]):
    stdrg += lopt.res.stripres.drag[i, 0]

trdrg = lopt.res.acs.dirx*lopt.res.trres.trfrctot

print(f'stdrg = {stdrg:g}')
print(f'trdrg = {trdrg:g}')

stcdi = stdrg/lopt.res.qfs/lsys.sref
trcdi = trdrg/lopt.res.qfs/lsys.sref

print(f'stcdi = {stcdi:g}')
print(f'trcdi = {trcdi:g}')

# %%
