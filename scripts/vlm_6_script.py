#%%
# Import Dependencies
from IPython.display import display
from matplotlib.pyplot import figure
from numpy import asarray, asin, atan2, degrees, pi, radians
from pygeom.geom1d import LinearSpline1D
from pyvlm import LatticeOptimum, LatticeResult, LatticeSystem
from pyvlm.tools import Bell, Mass, bell_lift_force_distribution
from pyvlm.tools.trim import LevelTrim

#%%
# Create Lattice System
jsonfilepath = '../files/Test_rhofw.json'
lsys = LatticeSystem.from_json(jsonfilepath)
display(lsys)

#%%
# Create Trim Scenario
m = 6.577089 # kg
mass = Mass('Test Mass', m, lsys.rref.x, lsys.rref.y, lsys.rref.z)
CL = 0.6
rho = 1.1448386410124347 # kg/m**3 at 704.3m altitude

trm = LevelTrim('CL = 0.6', lsys)
trm.set_density(rho)
trm.set_mass(mass)
trm.trim_speed_from_CL(CL)
display(trm)

#%%
#  Copy Lattice System
lsys_opt = lsys.copy_from_source()
lsys_bll = lsys.copy_from_source()

#%%
# Create Lattice Result
lres_org = trm.create_trim_result()
lres_org.name = 'Initial'
display(lres_org)

#%%
# Create Bell Lift Distribution
bll = Bell(lsys.bref, lsys.srfcs[0].strpy)
bll.set_density(trm.density)
bll.set_speed(trm.speed)
bll.set_lift(trm.lift)

lbll = bll.lift_force_distribution()

l = bell_lift_force_distribution(lsys.srfcs[0].strpy, lsys.bref, trm.lift)

#%%
# Bell Shaped Lift Distribution
lopt_bll = LatticeOptimum('Bell', lsys_bll)
lopt_bll.set_target_lift_force_distribution(l, trm.density, trm.speed)
lopt_bll.add_record('l', strplst=lsys.mstrpi)
display(lopt_bll)

#%%
# Phi Distribution Plot
axp = None
axp = lres_org.plot_phi_distribution(ax=axp)
axp = lopt_bll.plot_phi_distribution(ax=axp)
_ = axp.set_ylabel('Phi Distribution')
_ = axp.set_xlabel('Span Position')

#%%
# Lift Distribution Plot
axl = None
axl = lres_org.plot_trefftz_lift_force_distribution(ax=axl)
axl = lres_org.plot_strip_lift_force_distribution(ax=axl)
axl = lopt_bll.plot_trefftz_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%%
# Drag Distribution Plot
axd = None
axd = lres_org.plot_trefftz_drag_force_distribution(ax=axd)
axd = lres_org.plot_strip_drag_force_distribution(ax=axd)
axd = lopt_bll.plot_trefftz_drag_force_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

#%%
# Wash Distribution Plot
axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lopt_bll.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%%
# Optimal Strip Twist
al_bll = lopt_bll.optimum_strip_twist(crit=1e-2)

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
w = [trm.weight/rho/trm.speed/lsys_bll.bref*3/2*((2*yi/lsys_bll.bref)**2-0.5) for yi in yspec]

print(min(w)/min(lopt_bll.trres.trwsh))

#%%
# Reset Results
lres_org.reset()
lopt_bll.reset()

#%%
# Lift Distribution Plot
axl = None
axl = lres_org.plot_trefftz_lift_force_distribution(ax=axl)
axl = lres_org.plot_strip_lift_force_distribution(ax=axl)
axl = lopt_bll.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt_bll.plot_strip_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%%
# Drag Distribution Plot
axd = None
axd = lres_org.plot_trefftz_drag_force_distribution(ax=axd)
axd = lres_org.plot_strip_drag_force_distribution(ax=axd)
axd = lopt_bll.plot_trefftz_drag_force_distribution(ax=axd)
axd = lopt_bll.plot_strip_drag_force_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

#%%
# Wash Distribution Plot
axw = None
axw = lres_org.plot_trefftz_wash_distribution(ax=axw)
axw = lopt_bll.plot_trefftz_wash_distribution(ax=axw)
axw.plot(yspec, w, label='Bell Wash')
axw.legend()
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%%
# Geometruic Twist Distribution Plot
axa = lopt_bll.plot_strip_twist_distribution()
# axa.plot(lsys.srfcs[0].strpy, al_bll, label='alpha Bell')
axa.plot(yspec, alspec, label='alpha Specified')
leg = axa.legend()

#%%
# Induced Angle
alf = LinearSpline1D(asarray(yspec), asarray(alspec))
als = alf.evaluate_points_at_t(asarray(lsys.srfcs[0].strpy))

ali = [degrees(atan2(lopt_bll.trres.trwsh[i], trm.speed)) for i in range(len(lopt_bll.trres.trwsh))]
al0 = [1.0-abs(yi)*2/lsys_bll.bref for yi in lsys.srfcs[0].strpy]

alp = [alsj+alij+al0j for alsj, alij, al0j in zip(als, ali, al0)]

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(lsys.srfcs[0].strpy, ali, label='Induced Angle')
# ax.plot(lsys.srfcs[0].strpy, al0, label='Zero Lift Angle')
ax.plot(lsys.srfcs[0].strpy, als, label='Specified Angle')
ax.plot(lsys.srfcs[0].strpy, alp, label='Total Angle')
leg = ax.legend()

#%%
# Local Lift
cmax = 0.40005
cmin = 0.10008

c = [cmax-abs(yi)*2/lsys_bll.bref*(cmax-cmin) for yi in lsys.srfcs[0].strpy]

cla = 2*pi*0.91

l = [cla*ci*radians(alpi)*trm.dynpres for ci, alpi in zip(c, alp)]

fig  = figure(figsize=(12, 8))
ax = fig.gca()
ax.grid(True)
ax.plot(lsys.srfcs[0].strpy, l, label='Lift Distribution')
ax = lopt_bll.plot_trefftz_lift_force_distribution(ax=ax)
ax = lopt_bll.plot_strip_lift_force_distribution(ax=ax)
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
# New Results
lres_1g = LatticeResult('1g Result', lsys_bll)
lres_1g.set_state(alpha=0.0, speed=trm.speed)
lres_1g.set_density(rho=rho)
display(lres_1g)

CL0 = lres_1g.nfres.CL
CLa = lres_1g.stres.alpha.CL
al3g = degrees(asin((3*CL-CL0)/CLa))

lres_3g = LatticeResult('3g Result', lsys_bll)
lres_3g.set_state(alpha=al3g, speed=trm.speed)
lres_3g.set_density(rho=rho)
display(lres_3g)

#%%
# New Plots
axl = None
axl = lres_1g.plot_trefftz_lift_force_distribution(ax=axl)
axl = lres_1g.plot_strip_lift_force_distribution(ax=axl)
axl = lres_3g.plot_trefftz_lift_force_distribution(ax=axl)
axl = lres_3g.plot_strip_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

axd_1g = None
axd_1g = lres_1g.plot_trefftz_drag_force_distribution(ax=axd_1g)
axd_1g = lres_1g.plot_strip_drag_force_distribution(ax=axd_1g)
_ = axd_1g.set_ylabel('Drag Distribution')
_ = axd_1g.set_xlabel('Span Position')

axd_3g = None
axd_3g = lres_3g.plot_trefftz_drag_force_distribution(ax=axd_3g)
axd_3g = lres_3g.plot_strip_drag_force_distribution(ax=axd_3g)
_ = axd_3g.set_ylabel('Drag Distribution')
_ = axd_3g.set_xlabel('Span Position')

axw = None
axw = lres_1g.plot_trefftz_wash_distribution(ax=axw)
axw = lres_3g.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%%
# Display Strip Geometry
display(lsys_bll.strip_geometry)
