#%% Load Dependencies
from IPython.display import display
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm import latticesystem_from_json

#%% Create Lattice System
jsonfilepath = r'..\files\Straight_Wing_Cosine_100.json'
lsys = latticesystem_from_json(jsonfilepath)
display(lsys)

#%% Original Strip Geometry
display(lsys.strip_geometry)

#%% Original Strip Geometry
display(lsys.panel_geometry)

#%% Original Case

alpha = 3.0 # degrees
beta = 0.0 # degrees

lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=alpha, beta=beta)
display(lres_org)

# lres_org.print_total_loads()
# lres_org.print_aerodynamic_coefficients()

#%% Original Strip Forces
display(lres_org.strip_forces)

#%% Original Strip Coefficients
display(lres_org.strip_coefficients)

#%% Original Panel Forces
display(lres_org.panel_forces)

#%% Equivalent Elliptical Lift
from pyvlm.tools import elliptical_lift_distribution

L = lres_org.nfres.CL*lres_org.qfs*lsys.sref

lell = elliptical_lift_distribution(lsys.srfcs[0].strpy, lsys.bref, L)

lres_ell = LatticeResult('Equivalent Elliptical', lsys)
lres_ell.set_state(alpha=alpha)
lres_ell.set_lift_distribution(lell, rho=1.0, speed=1.0)
display(lres_ell)

# lres_ell.print_total_loads()
# lres_ell.print_aerodynamic_coefficients()

#%% Calculate Unconstrained Optimum

# Lspec = lres_org.CL*lres_org.qfs*lsys.sref

# phi1, lam1, Di1, L1, l1 = lsys.optimum_lift_distribution(Lspec)

# print(f'Di1 = {Di1}')
# print(f'L1 = {L1}')
# print(f'l1 = {l1}')

lopt_opt1 = LatticeOptimum('Unconstrained Optimised', lsys)
lopt_opt1.set_conditions()
lopt_opt1.add_constraint('L', L)
lopt_opt1.add_record('l', strplst='Mirrored')
phi1, lam1 = lopt_opt1.optimum_lift_distribution()
display(lopt_opt1)

lres_opt1 = LatticeResult('Unconstrained Optimised', lsys)
lres_opt1.set_state(alpha=alpha)
lres_opt1.set_phi(phi1)
display(lres_opt1)

# lres_opt1.print_total_loads()
# lres_opt1.print_aerodynamic_coefficients()

#%% Calculate Constrained Optimum

lspec = lopt_opt1.record[0].evaluate(lopt_opt1.pmat)*3/4

# phi2, lam2, Di2, L2, l2 = lsys.optimum_lift_distribution(L, lspec=lspec)

lopt_opt2 = LatticeOptimum('Constrained Optimised', lsys)
lopt_opt2.set_conditions()
lopt_opt2.add_constraint('L', L)
lopt_opt2.add_constraint('l', lspec, strplst='Mirrored')
phi2, lam2 = lopt_opt2.optimum_lift_distribution()
display(lopt_opt2)

# print(f'Di2 = {Di2}')
# print(f'L2 = {L2}')
# print(f'l2 = {l2}')

lres_opt2 = LatticeResult('Constrained Optimised', lsys)
lres_opt2.set_state(alpha=alpha)
lres_opt2.set_phi(phi2)
display(lres_opt2)

# lres_opt2.print_total_loads()
# lres_opt2.print_aerodynamic_coefficients()

#%% Plot Distribution

axl = None
axl = lres_org.plot_trefftz_lift_distribution(ax=axl)
axl = lres_ell.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt1.plot_trefftz_lift_distribution(ax=axl)
axl = lres_opt2.plot_trefftz_lift_distribution(ax=axl)

#%% Print Results

from math import pi

CDi_org_theory = lres_org.trres.CL**2/pi/lsys.ar/lres_org.trres.e

CDi_ell_theory = lres_org.trres.CL**2/pi/lsys.ar

print(f'CL_org = {lres_org.trres.CL:.3f}')
print(f'CL_ell = {lres_ell.trres.CL:.3f}')
print(f'CL_opt1 = {lres_opt1.trres.CL:.3f}')
print(f'CL_opt2 = {lres_opt2.trres.CL:.3f}')
print('')
print(f'CDi_org_theory = {CDi_org_theory:.7f}')
print(f'CDi_org = {lres_org.trres.CDi:.7f}')
print(f'CDi_ell_theory = {CDi_ell_theory:.7f}')
print(f'CDi_ell = {lres_ell.trres.CDi:.7f}')
print(f'CDi_opt1 = {lres_opt1.trres.CDi:.7f}')
print(f'CDi_opt2 = {lres_opt2.trres.CDi:.7f}')
print('')
print(f'Efficiency Improvement = {100.0*(1.0-lres_opt1.trres.CDi/lres_org.trres.CDi):.2f}%')
