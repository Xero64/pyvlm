#%% Load Dependencies
from pyvlm import latticesystem_from_json, LatticeSystem
from pyvlm.classes.latticesystem import velocity_matrix
from numpy.linalg import norm
from time import perf_counter

#%% Create Lattice System
jsonfilepath = r'..\files\Test_matrix.json'
lsys = latticesystem_from_json(jsonfilepath, build=False)
print(lsys)

#%% Calculate Matrices
_ = lsys.ra
_ = lsys.rb
_ = lsys.rc
_ = lsys.rg
_ = lsys.avc
_ = lsys.avg

#%% Define Functions
def control_velocity_matrix(lsys: LatticeSystem):
    from math import pi
    veli, vela, velb = velocity_matrix(lsys.ra, lsys.rb, lsys.rc)
    return (veli+vela-velb)/(4*pi)

def induced_velocity_matrix(lsys: LatticeSystem):
    from math import pi
    from numpy import fill_diagonal
    veli, vela, velb = velocity_matrix(lsys.ra, lsys.rb, lsys.rg)
    fill_diagonal(veli.x, 0.0)
    fill_diagonal(veli.y, 0.0)
    fill_diagonal(veli.z, 0.0)
    return (veli+vela-velb)/(4*pi)

#%% Build Control Velocity Matrix
start = perf_counter()
velc = control_velocity_matrix(lsys)
finish = perf_counter()
elapsed = finish-start
print(f'Built Control Velocity Matrix in {elapsed:.3f} seconds.')

#%% Check Control Difference
diffc = lsys.avc-velc
nrmcx = norm(diffc.x)
print(f'nrmc.x = {nrmcx}')
nrmcy = norm(diffc.y)
print(f'nrmc.y = {nrmcy}')
nrmcz = norm(diffc.z)
print(f'nrmc.z = {nrmcz}')

#%% Build Induced Velocity Matrix
start = perf_counter()
veli = induced_velocity_matrix(lsys)
finish = perf_counter()
elapsed = finish-start
print(f'Built Induced Velocity Matrix in {elapsed:.3f} seconds.')

#%% Check Induced Difference
diffi = lsys.avg-veli
nrmix = norm(diffi.x)
print(f'nrmi.x = {nrmix}')
nrmiy = norm(diffi.y)
print(f'nrmi.y = {nrmiy}')
nrmiz = norm(diffi.z)
print(f'nrmi.z = {nrmiz}')
