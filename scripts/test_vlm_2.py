#%% Import Dependencies

from pyvlm import latticesystem_from_json, LatticeResult, LatticeOptimum

#%% Low AR Wing

jsonfilepath1 = r"packages\pyvlm\files\Test_sweep.json"
lsys1 = latticesystem_from_json(jsonfilepath1)

print(f'bref 1 = {lsys1.bref:g}')
print(f'cref 1 = {lsys1.cref:g}')
print(f'sref 1 = {lsys1.sref:g}')
print(f'rref 1 = {lsys1.rref:.3f}')

#%% High AR Wing

jsonfilepath2 = r"packages\pyvlm\files\Test_sweep_2.json"
lsys2 = latticesystem_from_json(jsonfilepath2)

print(f'bref 2 = {lsys2.bref:g}')
print(f'cref 2 = {lsys2.cref:g}')
print(f'sref 2 = {lsys2.sref:g}')
print(f'rref 2 = {lsys2.rref:.3f}')

#%% Low AR Wing Optimum

lopt1 = LatticeOptimum('Low AR Wing', lsys1)
lopt1.set_conditions()
lopt1.add_constraint('L', 1.0)
lopt1.add_record('l', strplst='Mirrored')

phi1, lam1 = lopt1.optimum_lift_distribution()

lopt1.print_report()

lres1 = LatticeResult('Low AR Wing', lsys1)
lres1.set_conditions()
lres1.set_phi(phi1)
lres1.print_aerodynamic_coefficients()

#%% High AR Wing Constrained Optimum

l = lopt1.record[0].evaluate(lopt1.pmat)

lopt2 = LatticeOptimum('High AR Wing Constrained', lsys2)
lopt2.set_conditions()
lopt2.add_constraint('L', 1.0)
lopt2.add_constraint('l', l, strplst='Mirrored')

phi2, lam2 = lopt2.optimum_lift_distribution()

lopt2.print_report()

lres2 = LatticeResult('High AR Wing Constrained', lsys2)
lres2.set_conditions()
lres2.set_phi(phi2)
lres2.print_aerodynamic_coefficients()

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
