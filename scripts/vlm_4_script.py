#%%
# Import Dependencies
from math import atan2, degrees, radians, pi
from IPython.display import display
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm import latticesystem_from_json
from pyvlm.tools import bell_lift_force_distribution
from matplotlib.pyplot import figure
from pygeom.geom1d import LinearSpline

#%%
# Parameters
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

#%%
# Low AR Wing System
jsonfilepath = '../files/Test_rhofw.json'
lsys = latticesystem_from_json(jsonfilepath)
lsys_opt = latticesystem_from_json(jsonfilepath)
lsys_bll = latticesystem_from_json(jsonfilepath)
display(lsys)

#%%
# Create Initial Result
lres_org = LatticeResult('Initial', lsys)
lres_org.set_state(speed=V)
lres_org.set_density(rho=rho)
display(lres_org)

#%%
# Bell Shaped Lift Distribution
lbll = bell_lift_force_distribution(lsys.srfcs[0].strpy, lsys.bref, W)

lopt_bll = LatticeOptimum('Bell', lsys_bll)
lopt_bll.set_target_lift_force_distribution(lbll, rho, V)
lopt_bll.add_record('L')
lopt_bll.add_record('l', strplst=lsys.lstrpi)
lopt_bll.add_record('l', strplst=lsys.mstrpi)
display(lopt_bll)

#%%
# Optimal Lift Distribution
lopt = LatticeOptimum('Optimal', lsys_opt)
lopt.set_state(speed=V)
lopt.set_density(rho=rho)
lopt.add_constraint('L', W)
lopt.add_constraint('l', lopt_bll.record[1].value, strplst=lsys.lstrpi)
lopt.add_constraint('l', lopt_bll.record[2].value, strplst=lsys.mstrpi)
lopt.optimum_lift_force_distribution()
display(lopt)

#%%
# Plot Phi Distribution
axp = None
axp = lres_org.plot_phi_distribution(ax=axp)
axp = lopt_bll.plot_phi_distribution(ax=axp)
axp = lopt.plot_phi_distribution(ax=axp)
_ = axp.set_ylabel('Phi Distribution')
_ = axp.set_xlabel('Span Position')

#%%
# Plot Lift Distribution
axl = None
axl = lres_org.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt_bll.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt.plot_trefftz_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%%
# Plot Drag Distribution
axd = None
axd = lres_org.plot_trefftz_drag_force_distribution(ax=axd)
axd = lopt_bll.plot_trefftz_drag_force_distribution(ax=axd)
axd = lopt.plot_trefftz_drag_force_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

#%%
# Plot Wash Distribution
axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lopt_bll.plot_trefftz_wash_distribution(ax=axw)
axw = lopt.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%%
# Optimal Strip Twist
al_opt = lopt.optimum_strip_twist(crit=1e-1)
al_bll = lopt_bll.optimum_strip_twist(crit=1e-1)

#%%
# Specified String Twist

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

#%%
# Bell Downwash

w = [W/rho/V/lsys_bll.bref*3/2*((2*yi/lsys_bll.bref)**2-0.5) for yi in yspec]

print(min(w)/min(lopt_bll.trres.trwsh))

#%%
# Plot Lift Distribution
axl = None
axl = lres_org.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt_bll.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt.plot_strip_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%%
# Plot Drag Distribution
axd = None
axd = lres_org.plot_trefftz_drag_force_distribution(ax=axd)
axd = lopt_bll.plot_trefftz_drag_force_distribution(ax=axd)
axd = lopt.plot_trefftz_drag_force_distribution(ax=axd)
axd = lopt.plot_strip_drag_force_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

#%%
# Plot Wash Distribution
axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lopt_bll.plot_trefftz_wash_distribution(ax=axw)
axw = lopt.plot_trefftz_wash_distribution(ax=axw)
axw.plot(yspec, w, label='Bell Wash')
axw.legend()
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%%
# Plot Strip Geometric Twist
axa = lopt.plot_strip_twist_distribution()
axa.plot(lsys.srfcs[0].strpy, al_bll, label='alpha Bell')
axa.plot(lsys.srfcs[0].strpy, al_opt, label='alpha Optimum')
axa.plot(yspec, alspec, label='alpha Specified')
leg = axa.legend()

#%%
# Induced Angle
alf = LinearSpline(yspec, alspec)
als = alf.list_interpolation(lsys.srfcs[0].strpy)

ali = [degrees(atan2(lopt_bll.trres.trwsh[i], V)) for i in range(len(lopt_bll.trres.trwsh))]
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

#%%
# Local Lift

cmax = 0.40005
cmin = 0.10008

c = [cmax-abs(yi)*2/lsys_bll.bref*(cmax-cmin) for yi in lsys.srfcs[0].strpy]

cla = 2*pi*0.91

l = [cla*ci*radians(alpi)*q for ci, alpi in zip(c, alp)]

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(lsys.srfcs[0].strpy, l, label='Lift Distribution')
ax = lopt_bll.plot_trefftz_lift_force_distribution(ax=ax)
leg = ax.legend()

#%%
# Plot Specified Twist
fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(yspec, alspec)
ax.set_title('Wing Geometric Twist as a function of Span')
ax.set_xlabel('Span Coordinate - y - [m]')
_ = ax.set_ylabel('Geometric Twist Angle - $\\alpha_g$ - [deg]')

#%%
# Sum of Drag
stdrg = lopt.stripres.drag.sum()
trdrg = lopt.acs.dirx*lopt.trres.trfrctot

print(f'stdrg = {stdrg:g}')
print(f'trdrg = {trdrg:g}')

stcdi = stdrg/lopt.qfs/lsys.sref
trcdi = trdrg/lopt.qfs/lsys.sref

print(f'stcdi = {stcdi:g}')
print(f'trcdi = {trcdi:g}')
