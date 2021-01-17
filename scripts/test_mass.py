#%% Import Dependencies
from pyvlm.tools import masses_from_json, MassCollection

#%% Read in Mass File
jsonfilepath = '../files/Aircraft_Mass.json'
masses = masses_from_json(jsonfilepath)
masscol = MassCollection('Aircraft Mass', masses)
print(masscol)
