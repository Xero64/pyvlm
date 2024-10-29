#%%
# Import Dependencies
from sympy import Symbol, cos, sin
from pygeom.symbol3d import SymbolicVector

#%%
# Create Symbols
theta = Symbol('theta', real=True)

v = SymbolicVector(1, 0, 0)
k = SymbolicVector(0, 0, 1)

vrot: SymbolicVector = v*cos(theta) + k.cross(v)*sin(theta) + k.dot(v)*(1 - cos(theta))*k
vrot = vrot.simplify()
print(f'{vrot = }\n')

dvrot = vrot.diff(theta)
print(f'{dvrot = }\n')

Dvrot = vrot - v
Dvrot = Dvrot.simplify()
print(f'{Dvrot = }\n')
