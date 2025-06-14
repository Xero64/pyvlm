'''
This script performs the span load optimization.
The low aspect ratio wing gives an elliptical spanload under optimisation.
The high aspect ratio wing is constrained to give the same root bending moment,
as that for the low aspect ratio wing resulting in the R.T. Jones distribution.
'''

#%%
# Import Dependencies
from IPython.display import display
from pyvlm import LatticeOptimum, LatticeSystem

#%%
# Low AR Wing
jsonfilepath1 = '../files/Sweep_Low_AR_100.json'
lsys1 = LatticeSystem.from_json(jsonfilepath1)
display(lsys1)

#%%
# High AR Wing
jsonfilepath2 = '../files/Sweep_High_AR_100.json'
lsys2 = LatticeSystem.from_json(jsonfilepath2)
display(lsys2)

#%%
# Low AR Wing Optimum
lopt1 = LatticeOptimum('Low AR Wing', lsys1)
lopt1.set_state()
lopt1.add_constraint('L', 1.0)
lopt1.add_record('l', strplst=lsys1.lstrpi)
lopt1.add_record('l', strplst=lsys1.mstrpi)
phi1, lam1 = lopt1.optimum_lift_force_distribution()
display(lopt1)

#%%
# High AR Wing Constrained Optimum
lopt2 = LatticeOptimum('High AR Wing Constrained', lsys2)
lopt2.set_state()
lopt2.add_constraint('L', 1.0)
lopt2.add_constraint('l', lopt1.record[0].value, strplst=lsys2.lstrpi)
lopt2.add_constraint('l', lopt1.record[1].value, strplst=lsys2.mstrpi)
phi2, lam2 = lopt2.optimum_lift_force_distribution()
display(lopt2)

#%%
# Print Drag Ratio

Di1 = lopt1.return_induced_drag()
Di2 = lopt2.return_induced_drag()
print(f'Drag Ratio = {Di2/Di1*100:.2f}%')

#%%
# Lift Distribution Plot
axl = None
axl = lopt1.plot_trefftz_lift_force_distribution(ax=axl)
axl = lopt2.plot_trefftz_lift_force_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

#%%
# Drag Distribution Plot
axd = None
axd = lopt1.plot_trefftz_drag_force_distribution(ax=axd)
axd = lopt2.plot_trefftz_drag_force_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

#%%
# Wash Distribution Plot
axw = None
axw = lopt1.plot_trefftz_wash_distribution(ax=axw)
axw = lopt2.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')
