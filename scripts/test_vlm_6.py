#%% Import Dependencies
from IPython.display import display_markdown
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm import latticesystem_from_json
from pyvlm.tools import bell_lift_distribution
from pyvlm.tools.trim import LevelTrim
from pyvlm.tools import Bell

#%% Create Lattice System

jsonfilepath = r'..\files\Test_rhofw.json'
lsys = latticesystem_from_json(jsonfilepath)
lsys.build()
display_markdown(lsys)

#%% Create Trim Scenario
m = 6.577089 # kg
CL = 0.6
rho = 1.1448386410124347 # kg/m**3 at 704.3m altitude

trm = LevelTrim('CL = 0.6', lsys)
trm.set_density(rho)
trm.set_mass(m)
trm.trim_speed_from_CL(CL)
display_markdown(trm)

#%%  Copy Lattice System
lsys_opt = lsys.copy_from_source()
lsys_bll = lsys.copy_from_source()
# lsys_opt = latticesystem_from_json(jsonfilepath)
# lsys_bll = latticesystem_from_json(jsonfilepath)

#%% Create Lattice Result
lres_org = trm.create_trim_result()
lres_org.name = 'Initial'
display_markdown(lres_org)

#%% Create Bell Lift Distribution
bll = Bell(lsys.bref, lsys.srfcs[0].strpy)
bll.set_density(trm.density)
bll.set_speed(trm.speed)
bll.set_lift(trm.lift)

lbll = bll.lift_distribution()

l = bell_lift_distribution(lsys.srfcs[0].strpy, lsys.bref, trm.lift)

#%% Bell Shaped Lift Distribution

lopt_bll = LatticeOptimum('Bell', lsys_bll)
lopt_bll.set_lift_distribution(l, trm.density, trm.speed)
lopt_bll.add_record('l', strplst='Mirrored')
display_markdown(lopt_bll)

lres_bll = LatticeResult('Bell', lsys)
lres_bll.set_lift_distribution(l, trm.density, trm.speed)
display_markdown(lres_bll)

#%% Plots

axp = None
axp = lres_org.plot_phi_distribution(ax=axp)
axp = lres_bll.plot_phi_distribution(ax=axp)
_ = axp.set_ylabel('Phi Distribution')
_ = axp.set_xlabel('Span Position')

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_org.plot_strip_lift_distribution(ax=axl)
axl = lres_bll.plot_trefftz_lift_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lres_org.plot_strip_drag_distribution(ax=axd)
axd = lres_bll.plot_trefftz_drag_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lres_bll.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%% Optimal Strip Twist

al_bll = lopt_bll.optimum_strip_twist(crit=1e-2)

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

# w = bll.trefftz_wash_distribution()

w = [trm.weight/rho/trm.speed/lsys_bll.bref*3/2*((2*yi/lsys_bll.bref)**2-0.5) for yi in yspec]

print(min(w)/min(lres_bll.trres.trwsh))

#%% Plots

lres_org.reset()
lres_bll.reset()

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_org.plot_strip_lift_distribution(ax=axl)
axl = lres_bll.plot_trefftz_lift_distribution(ax=axl)
axl = lres_bll.plot_strip_lift_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

axd = None
axd = lres_org.plot_trefftz_drag_distribution(ax=axd)
axd = lres_org.plot_strip_drag_distribution(ax=axd)
axd = lres_bll.plot_trefftz_drag_distribution(ax=axd)
axd = lres_bll.plot_strip_drag_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lres_bll.plot_trefftz_wash_distribution(ax=axw)
axw.plot(yspec, w, label='Bell Wash')
axw.legend()
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

axa = lopt_bll.plot_strip_twist_distribution()
# axa.plot(lsys.srfcs[0].strpy, al_bll, label='alpha Bell')
axa.plot(yspec, alspec, label='alpha Specified')
leg = axa.legend()

#%% Induced Angle

from math import atan2, degrees, radians, pi, asin
from matplotlib.pyplot import figure
from pygeom.geom1d import LinearSpline

alf = LinearSpline(yspec, alspec)
als = alf.list_interpolation(lsys.srfcs[0].strpy)

ali = [degrees(atan2(lres_bll.trres.trwsh[i], trm.speed)) for i in range(len(lres_bll.trres.trwsh))]
# al0 = [1.0-abs(yi)*2/lsys_bll.bref for yi in lsys.srfcs[0].strpy]

alp = [als[i]+ali[i] for i in range(len(lsys.srfcs[0].strpy))]

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(lsys.srfcs[0].strpy, ali, label='Induced Angle')
# ax.plot(lsys.srfcs[0].strpy, al0, label='Zero Lift Angle')
ax.plot(lsys.srfcs[0].strpy, als, label='Specified Angle')
ax.plot(lsys.srfcs[0].strpy, alp, label='Total Angle')
leg = ax.legend()

#%% Local Lift

cmax = 0.40005
cmin = 0.10008

c = [cmax-abs(yi)*2/lsys_bll.bref*(cmax-cmin) for yi in lsys.srfcs[0].strpy]

cla = 2*pi*0.91

l = [cla*c[i]*radians(alp[i])*trm.dynpres for i in range(len(lsys.srfcs[0].strpy))]

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(lsys.srfcs[0].strpy, l, label='Lift Distribution')
ax = lres_bll.plot_trefftz_lift_distribution(ax=ax)
ax = lres_bll.plot_strip_lift_distribution(ax=ax)
leg = ax.legend()

#%% Plot Specified Twist

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(yspec, alspec)
ax.set_title('Wing Geometric Twist as a function of Span')
ax.set_xlabel('Span Coordinate - y - [m]')
_ = ax.set_ylabel('Geometric Twist Angle - $\\alpha_g$ - [deg]')

#%% New Results

lres_1g = LatticeResult('1g Result', lsys_bll)
lres_1g.set_state(alpha=0.0, speed=trm.speed)
lres_1g.set_density(rho=rho)
display_markdown(lres_1g)

CL0 = lres_1g.nfres.CL
CLa = lres_1g.stres['alpha'].CL
al3g = degrees(asin((3*CL-CL0)/CLa))

lres_3g = LatticeResult(f'3g Result', lsys_bll)
lres_3g.set_state(alpha=al3g, speed=trm.speed)
lres_3g.set_density(rho=rho)
display_markdown(lres_3g)

#%% New Plots

axl = None
axl = lres_1g.plot_trefftz_lift_distribution(ax=axl)
axl = lres_1g.plot_strip_lift_distribution(ax=axl)
axl = lres_3g.plot_trefftz_lift_distribution(ax=axl)
axl = lres_3g.plot_strip_lift_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

axd0 = None
axd0 = lres_1g.plot_trefftz_drag_distribution(ax=axd0)
axd0 = lres_1g.plot_strip_drag_distribution(ax=axd0)
_ = axd0.set_ylabel('Drag Distribution')
_ = axd0.set_xlabel('Span Position')

axd14 = None
axd14 = lres_3g.plot_trefftz_drag_distribution(ax=axd14)
axd14 = lres_3g.plot_strip_drag_distribution(ax=axd14)
_ = axd14.set_ylabel('Drag Distribution')
_ = axd14.set_xlabel('Span Position')

axw = None
axw = lres_1g.plot_trefftz_wash_distribution(ax=axw)
axw = lres_3g.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%% Display Strip Geometry
display_markdown(lsys_bll.strip_geometry)
