#%% Import Dependencies
from pyvlm.tools import masses_from_json
from pyvlm.tools.mass import display_masses

#%% Read in Mass File
jsonfilepath = r'..\files\Aircraft_Mass.json'
masses = masses_from_json(jsonfilepath)

#%% Display Masses
display_masses(masses)
