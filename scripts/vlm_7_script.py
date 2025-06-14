#%%
# Load Dependencies
from IPython.display import display
from pyvlm import LatticeOptimum, LatticeSystem
from pyvlm.tools import elliptical_lift_force_distribution

#%%
# Create Lattice System
jsonfilepath = '../files/Straight_Wing_Cosine_b100_c1.json'
lsys = LatticeSystem.from_json(jsonfilepath)
display(lsys)

#%%
# Create Elliptical Distribution
lell = elliptical_lift_force_distribution(lsys.srfcs[0].strpy, lsys.bref, 10000.0)

#%%
# Create Elliptical Optimum
lopt_ell = LatticeOptimum('Equivalent Elliptical', lsys)
lopt_ell.set_target_lift_force_distribution(lell, rho=1.0, speed=100.0)
display(lopt_ell)

#%%
# Determine Optimal Twist
al_ell = lopt_ell.optimum_strip_twist(crit=1e-6)

#%%
# Plot Optimum Twist
axs = lopt_ell.plot_strip_twist_distribution()

#%%
# Plot Optimal Wash
axw = None
axw = lopt_ell.plot_trefftz_wash_distribution(ax=axw)
axw = lopt_ell.plot_panel_near_field_velocities(ax=axw, component='z')

#%%
# Print Near Field Results
display(lopt_ell.panel_near_field_results)
