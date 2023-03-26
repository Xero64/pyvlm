#%%
# Load Dependencies
from math import pi
from IPython.display import display
from pyvlm import LatticeResult, LatticeOptimum
from pyvlm import latticesystem_from_json
from pyvlm.tools import elliptical_lift_force_distribution

#%%
# Create Lattice System
jsonfilepath = '../files/Straight_Wing_Cosine_100.json'
lsys = latticesystem_from_json(jsonfilepath)
display(lsys)

#%%
# Original Case
alpha = 3.0 # degrees
beta = 0.0 # degrees

lres_org = LatticeResult('Baseline', lsys)
lres_org.set_state(alpha=alpha, beta=beta)
display(lres_org)

#%%
# Equivalent Elliptical Lift
L = lres_org.nfres.CL*lres_org.qfs*lsys.sref
lell = elliptical_lift_force_distribution(lsys.srfcs[0].strpy, lsys.bref, L)

lopt_ell = LatticeOptimum('Equivalent Elliptical', lsys)
lopt_ell.set_state(alpha=alpha)
lopt_ell.set_target_lift_force_distribution(lell, rho=1.0, speed=1.0)
display(lopt_ell)

#%%
# Calculate Unconstrained Optimum
# Constrain Only Lift and the resulting lift distribution is Elliptical.
# Note: You have to consider root moment for both sides of the lifting surface.

lopt1 = LatticeOptimum('Unconstrained Optimised', lsys)#, sym=False)
lopt1.set_state()
lopt1.add_constraint('L', L)
lopt1.add_record('l', strplst=lsys.lstrpi)
lopt1.add_record('l', strplst=lsys.mstrpi)
phi1, lam1 = lopt1.optimum_lift_force_distribution()
display(lopt1)

#%%
# Calculate Constrained Optimum
# Constrain Root Bending Moment (Rolling Moment) to 75% of the Elliptical Optimimum.
# Note: You have to consider both sides of the lifting surface.

lspec = lopt1.record[0].value*0.75
mspec = lopt1.record[1].value*0.75

lopt2 = LatticeOptimum('Constrained Optimised', lsys)#, sym=False)
lopt2.set_state()
lopt2.add_constraint('L', L)
lopt2.add_constraint('l', lspec, strplst=lsys.lstrpi)
lopt2.add_constraint('l', mspec, strplst=lsys.mstrpi)
phi2, lam2 = lopt2.optimum_lift_force_distribution()
display(lopt2)

#%%
# Plot Distributions
axl = None
axl = lres_org.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt_ell.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt1.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt2.plot_trefftz_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%%
# Print Results
CDi_org_theory = lres_org.trres.CL**2/pi/lsys.ar/lres_org.trres.e
CDi_ell_theory = lres_org.trres.CL**2/pi/lsys.ar

print(f'CL_org = {lres_org.trres.CL:.3f}')
print(f'CL_ell = {lopt_ell.trres.CL:.3f}')
print(f'CL_opt1 = {lopt1.trres.CL:.3f}')
print(f'CL_opt2 = {lopt2.trres.CL:.3f}')
print('')
print(f'CDi_org_theory = {CDi_org_theory:.7f}')
print(f'CDi_org = {lres_org.trres.CDi:.7f}')
print(f'CDi_ell_theory = {CDi_ell_theory:.7f}')
print(f'CDi_ell = {lopt_ell.trres.CDi:.7f}')
print(f'CDi_opt1 = {lopt1.trres.CDi:.7f}')
print(f'CDi_opt2 = {lopt2.trres.CDi:.7f}')
print('')
print(f'Efficiency Improvement = {100.0*(1.0-lopt1.trres.CDi/lres_org.trres.CDi):.2f}%')
