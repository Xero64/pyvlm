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
display_markdown(lsys)

#%% Design Point Result
lres = lsys.results['Design Point']
display_markdown(lres)

#%% Plots
axp = None
axp = lres.plot_phi_distribution(ax=axp)
_ = axp.set_ylabel('Phi Distribution')
_ = axp.set_xlabel('Span Position')

axl = None
axl = lres.plot_trefftz_lift_distribution(ax=axl)
axl = lres.plot_strip_lift_distribution(ax=axl)
_ = axl.set_ylabel('Lift Distribution')
_ = axl.set_xlabel('Span Position')

axd = None
axd = lres.plot_trefftz_drag_distribution(ax=axd)
axd = lres.plot_strip_drag_distribution(ax=axd)
_ = axd.set_ylabel('Drag Distribution')
_ = axd.set_xlabel('Span Position')

axw = None
axw = lres.plot_trefftz_wash_distribution(ax=axw)
_ = axw.set_ylabel('Wash Distribution')
_ = axw.set_xlabel('Span Position')

#%% Display Strip Geometry
display_markdown(lsys.strip_geometry)
