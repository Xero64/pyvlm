#%% Import Dependencies
from IPython.display import display_markdown
from pygeom.geom3d import Point
from pyvlm.classes import LatticeSystem
from pyvlm.tools.trim import TurningTrim

#%% Create Turning Trim
lsys = LatticeSystem('Test Trim', [], 20.0, 1.0, 20.0, Point(0.0, 0.0, 0.0))
display_markdown(lsys)
trim = TurningTrim('30deg Banked Turn', lsys)
trim.set_speed_and_density(25.0, 1.2)
trim.set_mass(80.0)
trim.set_bank_angle(30.0)
display_markdown(trim)
