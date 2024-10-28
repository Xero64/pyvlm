#%%
# Import Dependencies
from sympy import Symbol, cos, sin
from pygeom.symbol3d import SymbolicVector

#%%
# Create Symbols
alpha = Symbol('alpha', real=True)
beta = Symbol('beta', real=True)

#%%
# Create Vectors for Aerodynamic Coordinate System
dirx = SymbolicVector(cos(alpha)*cos(beta), -sin(beta), sin(alpha)*cos(beta))
diry = SymbolicVector(cos(alpha)*sin(beta), cos(beta), sin(alpha)*sin(beta))
dirz = SymbolicVector(-sin(alpha), 0, cos(alpha))
print(f'{dirx = }\n')
print(f'{diry = }\n')
print(f'{dirz = }\n')

dirz = dirx.cross(diry).simplify()
print(f'Check: {dirz = }\n')

ddirx_dalpha = dirx.diff(alpha)
ddiry_dalpha = diry.diff(alpha)
ddirz_dalpha = dirz.diff(alpha)
print(f'{ddirx_dalpha = }\n')
print(f'{ddiry_dalpha = }\n')
print(f'{ddirz_dalpha = }\n')

ddirz_dalpha = (ddirx_dalpha.cross(diry) + dirx.cross(ddiry_dalpha)).simplify()
print(f'Check: {ddirz_dalpha = }\n')

ddirx_dbeta = dirx.diff(beta)
ddiry_dbeta = diry.diff(beta)
ddirz_dbeta = dirz.diff(beta)
print(f'{ddirx_dbeta = }\n')
print(f'{ddiry_dbeta = }\n')
print(f'{ddirz_dbeta = }\n')

ddirz_dbeta = (ddirx_dbeta.cross(diry) + dirx.cross(ddiry_dbeta)).simplify()
print(f'Check: {ddirz_dbeta = }\n')

#%%
# Create Vectors for Stability Coordinate System
dirx = SymbolicVector(-cos(alpha), 0, -sin(alpha))
diry = SymbolicVector(0, 1, 0)
dirz = SymbolicVector(sin(alpha), 0, -cos(alpha))
print(f'{dirx = }\n')
print(f'{diry = }\n')
print(f'{dirz = }\n')

dirz = dirx.cross(diry).simplify()
print(f'Check: {dirz = }\n')

ddirx_dalpha = dirx.diff(alpha)
ddiry_dalpha = diry.diff(alpha)
ddirz_dalpha = dirz.diff(alpha)

print(f'{ddirx_dalpha = }\n')
print(f'{ddiry_dalpha = }\n')
print(f'{ddirz_dalpha = }\n')

ddirz_dalpha = (ddirx_dalpha.cross(diry) + dirx.cross(ddiry_dalpha)).simplify()
print(f'Check: {ddirz_dalpha = }\n')

ddirx_dbeta = dirx.diff(beta)
ddiry_dbeta = diry.diff(beta)
ddirz_dbeta = dirz.diff(beta)
print(f'{ddirx_dbeta = }\n')
print(f'{ddiry_dbeta = }\n')
print(f'{ddirz_dbeta = }\n')

ddirz_dbeta = (ddirx_dbeta.cross(diry) + dirx.cross(ddiry_dbeta)).simplify()
print(f'Check: {ddirz_dbeta = }\n')
