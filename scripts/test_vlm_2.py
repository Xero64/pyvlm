#%% Import Dependencies
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm.files import load_package_file

#%% Low AR Wing
jsonfilename1 = "Sweep_Low_AR_100.json"
lsys1 = load_package_file(jsonfilename1)
print(lsys1)

#%% High AR Wing
jsonfilename2 = "Sweep_High_AR_100.json"
lsys2 = load_package_file(jsonfilename2)
print(lsys2)

#%% Low AR Wing Optimum
lopt1 = LatticeOptimum('Low AR Wing', lsys1)
lopt1.set_conditions()
lopt1.add_constraint('L', 1.0)
lopt1.add_record('l', strplst='Mirrored')
phi1, lam1 = lopt1.optimum_lift_distribution()
print(lopt1)

lres1 = LatticeResult('Low AR Wing', lsys1)
lres1.set_conditions()
lres1.set_phi(phi1)
print(lres1)

#%% Constrained Root Bending Moment
l = lopt1.record[0].evaluate(lopt1.pmat)

#%% High AR Wing Constrained Optimum
lopt2 = LatticeOptimum('High AR Wing Constrained', lsys2)
lopt2.set_conditions()
lopt2.add_constraint('L', 1.0)
lopt2.add_constraint('l', l, strplst='Mirrored')
phi2, lam2 = lopt2.optimum_lift_distribution()
print(lopt2)

lres2 = LatticeResult('High AR Wing Constrained', lsys2)
lres2.set_conditions()
lres2.set_phi(phi2)
print(lres2)

#%% Print Drag Ratio

Di1 = lopt1.return_induced_drag()
Di2 = lopt2.return_induced_drag()
print(f'\nDrag Ratio = {Di2/Di1*100:.2f}%')

#%% Plots

axl = None
axl = lres1.plot_trefftz_lift_distribution(ax=axl)
axl = lres2.plot_trefftz_lift_distribution(ax=axl)

axd = None
axd = lres1.plot_trefftz_drag_distribution(ax=axd)
axd = lres2.plot_trefftz_drag_distribution(ax=axd)

axw = None
axw = lres1.plot_trefftz_wash_distribution(ax=axw)
axw = lres2.plot_trefftz_wash_distribution(ax=axw)
