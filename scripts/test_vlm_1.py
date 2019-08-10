#%% Load lattice JSON File
from pyvlm import LatticeResult
from pyvlm.files import load_package_file

jsonfilename = "Straight_Wing_Cosine_20.json"
lsys = load_package_file(jsonfilename)

print(f'bref = {lsys.bref:g}')
print(f'cref = {lsys.cref:g}')
print(f'sref = {lsys.sref:g}')
print(f'rref = {lsys.rref:.3f}')

#%% Original Case

alpha = 3.0 # degrees

lres_org = LatticeResult('Baseline', lsys)
lres_org.set_conditions(alpha=alpha)

# lres_org.print_total_loads()
lres_org.print_aerodynamic_coefficients()

#%% Original Strip Forces
lres_org.print_strip_forces()

#%% Equivalent Elliptical Lift
from pyvlm.tools import elliptical_lift_distribution

L = lres_org.CL*lres_org.qfs*lsys.sref

lell = elliptical_lift_distribution(lres_org.strpy, lsys.bref, L)

lres_ell = LatticeResult('Equivalent Elliptical', lsys)
lres_ell.set_conditions(alpha=alpha)
lres_ell.set_lift_distribution(lell, rho=1.0, speed=1.0)

# lres_ell.print_total_loads()
lres_ell.print_aerodynamic_coefficients()

#%% Calculate Unconstrained Optimum

Lspec = lres_org.CL*lres_org.qfs*lsys.sref

phi1, lam1, Di1, L1, l1 = lsys.optimum_lift_distribution(Lspec)

print(f'Di1 = {Di1}')
print(f'L1 = {L1}')
print(f'l1 = {l1}')

lres_opt1 = LatticeResult('Unconstrained Optimised', lsys)
lres_opt1.set_conditions(alpha=alpha)
lres_opt1.set_phi(phi1)

# lres_opt1.print_total_loads()
lres_opt1.print_aerodynamic_coefficients()

#%% Calculate Constrained Optimum

phi2, lam2, Di2, L2, l2 = lsys.optimum_lift_distribution(Lspec, lspec=l1*3/4)

print(f'Di2 = {Di2}')
print(f'L2 = {L2}')
print(f'l2 = {l2}')

lres_opt2 = LatticeResult('Constrained Optimised', lsys)
lres_opt2.set_conditions(alpha=alpha)
lres_opt2.set_phi(phi2)

# lres_opt2.print_total_loads()
lres_opt2.print_aerodynamic_coefficients()

#%% Plot Distribution

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_ell.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt1.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt2.plot_trefftz_lift_distribution(ax=axl)

#%% Print Results

from math import pi

CDi_org_theory = lres_org.CL_ff**2/pi/lsys.ar/lres_org.e

CDi_ell_theory = lres_org.CL_ff**2/pi/lsys.ar

print(f'CL_org = {lres_org.CL_ff:.7f}')
print(f'CL_ell = {lres_ell.CL_ff:.7f}')
print(f'CL_opt1 = {lres_opt1.CL_ff:.7f}')
print(f'CL_opt2 = {lres_opt2.CL_ff:.7f}')
print('')
print(f'CDi_org_theory = {CDi_org_theory:.7f}')
print(f'CDi_org = {lres_org.CDi_ff:.7f}')
print(f'CDi_ell_theory = {CDi_ell_theory:.7f}')
print(f'CDi_ell = {lres_ell.CDi_ff:.7f}')
print(f'CDi_opt1 = {lres_opt1.CDi_ff:.7f}')
print(f'CDi_opt2 = {lres_opt2.CDi_ff:.7f}')
print('')
print(f'Efficiency Improvement = {100.0*(1.0-lres_opt1.CDi_ff/lres_org.CDi_ff):.2f}%')

#%% Plot Surfaces

lsys.plot_surface(view='top')

# #%% Pause

# y = [pnt.y for pnt in lsys.srfcs[0].pnts[:, 0].transpose().tolist()[0]]

# from matplotlib.pyplot import show
# show()
# pass
