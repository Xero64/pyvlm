#%%
# Load Dependencies
from math import pi
from time import perf_counter
from pyvlm import latticesystem_from_json, LatticeSystem
from pyvlm.classes.latticesystem import velocity_matrix
from numpy.linalg import norm
from numpy import fill_diagonal

#%%
# Create Lattice System
jsonfilepath = '../files/Test_matrix.json'
lsys = latticesystem_from_json(jsonfilepath, mesh=False)
print(lsys)

#%%
# Calculate Matrices
_ = lsys.ra
_ = lsys.rb
_ = lsys.rc
_ = lsys.rg
_ = lsys.avc
_ = lsys.avg

#%%
# Define Functions
def control_velocity_matrix(sys: LatticeSystem):
    vgi, vga, vgb = velocity_matrix(sys.ra, sys.rb, sys.rc)
    return (vgi + vga - vgb)/(4*pi)

def induced_velocity_matrix(sys: LatticeSystem):
    vgi, vga, vgb = velocity_matrix(sys.ra, sys.rb, sys.rg,)
    fill_diagonal(vgi.x, 0.0)
    fill_diagonal(vgi.y, 0.0)
    fill_diagonal(vgi.z, 0.0)
    return (vgi + vga - vgb)/(4*pi)

#%%
# Build Induced Velocity Matrix
start = perf_counter()
veli = induced_velocity_matrix(lsys)
finish = perf_counter()
elapsed = finish-start
print(f'Built Induced Velocity Matrix in {elapsed:.3f} seconds.')

#%%
# Check Induced Difference
diffi = lsys.avg(0.0)-veli
nrmix = norm(diffi.x)
print(f'nrmi.x = {nrmix}')
nrmiy = norm(diffi.y)
print(f'nrmi.y = {nrmiy}')
nrmiz = norm(diffi.z)
print(f'nrmi.z = {nrmiz}')

#%%
# Build Control Velocity Matrix
start = perf_counter()
velc = control_velocity_matrix(lsys)
finish = perf_counter()
elapsed = finish-start
print(f'Built Control Velocity Matrix in {elapsed:.3f} seconds.')

#%%
# Check Control Difference
diffc = lsys.avc(0.0)-velc
nrmcx = norm(diffc.x)
print(f'nrmc.x = {nrmcx}')
nrmcy = norm(diffc.y)
print(f'nrmc.y = {nrmcy}')
nrmcz = norm(diffc.z)
print(f'nrmc.z = {nrmcz}')
