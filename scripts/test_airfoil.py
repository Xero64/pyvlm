#%% Import Dependencies
from pyvlm.tools import airfoil_from_dat, Airfoil
from matplotlib.pyplot import figure
from numpy.linalg import lstsq, solve
from numpy.matlib import zeros
from math import acos, floor, cos, sin, pi, sqrt
from pyfoil.airfoil import naca_to_xyt

#%% Create Airfoil
# datfilepath = r'C:\SSD\RhoFW\root_airfoil.dat'
# airfoil = airfoil_from_dat(datfilepath)

code = '0012'
x, y, t = naca_to_xyt(code)
airfoil = Airfoil(code, x, y)

# print(airfoil.rngt)
# print(airfoil.rngb)

# indt = [i for i in airfoil.rngt]
# indb = [i for i in airfoil.rngb]

# xnt = [airfoil.xn[i] for i in airfoil.rngt]
# xnb = [airfoil.xn[i] for i in airfoil.rngb]

# x = [1.0, 0.5, 0.0, 0.5, 1.0]
# y = [0.0, 0.1, 0.0, 0.1, 0.0]

#%% Plot Profile
ax1 = airfoil.plot_airfoil()
ax2 = airfoil.plot_normalised_aifoil()

#%% Fitting Kernels

def camber_kernel(n: int, th: float):
    if n == 1:
        return sin(th)**2/4
    else:
        return (n*sin(th)*sin(n*th) + cos(th)*cos(n*th) - 1)/(2*(n**2 - 1))

def camber_slope_kernel(n: int, th: float):
    return cos(n*th)

def thickness_kernel(n: int, th: float):
    return sin(n*thi)

def thickness_slope_kernel(n: int, th: float):
    return 2*n*cos(n*th)/sin(th)

# def camber_kernel(n: int, x: float):
#     return x**n

# def thickness_kernel(n: int, x: float):
#     if n == 0:
#         return sqrt(x)
#     else:
#         return x**n

# x = c/2*(1-cos(th))
# th = acos(1-2*x/c)
num = 12
xu, yu, xl, yl = airfoil.cosine_spline(num=num)

thu = [acos(1-2*xi) for xi in xu]
thl = [acos(1-2*xi) for xi in xl]

# tht = [acos(1.0-2.0*airfoil.xn[i]) for i in airfoil.rngt]
# thb = [acos(1.0-2.0*airfoil.xn[i]) for i in airfoil.rngb]

# num = 2*len(thb)+2*len(tht)
# num = len(xnt)+len(xnb)
# numt = len(thb)+len(tht)
# numc = len(thb)+len(tht)
# num = numt+numc

numu = len(thu)
numl = len(thl)
numc = int((numu+numl+2)/2)
# numt = int((numu+numl+2)/2)
numt = 0

# numl = 0
# numt = numt-numl
# numc = numc-numl

a = zeros((numu+numl+2, numc+numt))
b = zeros((numu+numl+2, 1))

#%% Fill Matrix Terms
i = 0
# Top Point Fit
for c, thi in enumerate(thu):
    j = 0
    for k in range(numc):
        n = k+1
        a[i, j] = camber_kernel(n, thi)
        j += 1
    for k in range(numt):
        n = k+1
        a[i, j] = thickness_kernel(n, thi)/2
        j += 1
    # ind = indt[c]
    b[i, 0] = yu[c]
    i += 1

# Bottom Point Fit
for c, thi in enumerate(thl):
    j = 0
    for k in range(numc):
        n = k+1
        a[i, j] = camber_kernel(n, thi)
        j += 1
    for k in range(numt):
        n = k+1
        a[i, j] = -thickness_kernel(n, thi)/2
        j += 1
    # ind = indb[c]
    b[i, 0] = yl[c]
    i += 1

j = 0
for k in range(numc):
    n = k+1
    a[i, j] = camber_kernel(n, pi)
    j += 1
b[i, 0] = airfoil.yte
i += 1

# Trailing Edge Camber Slope Fit
j = 0
for k in range(numc):
    n = k+1
    a[i, j] = camber_slope_kernel(n, pi)
    j += 1
b[i, 0] = airfoil.mte
i += 1

# Top Slope Fit
# for c, thi in enumerate(tht):
#     j = 0
#     for k in range(numc):
#         n = k+1
#         a[i, j] = camber_slope_kernel(n, thi)
#         j += 1
#     for k in range(numt):
#         n = k+1
#         a[i, j] = thickness_slope_kernel(n, thi)/2
#         j += 1
#     ind = indt[c]
#     b[i, 0] = airfoil.dydx[ind]
#     i += 1

# Bottom Slope Fit
# for c, thi in enumerate(thb):
#     j = 0
#     for k in range(numc):
#         n = k+1
#         a[i, j] = camber_slope_kernel(n, thi)
#         j += 1
#     for k in range(numt):
#         n = k+1
#         a[i, j] = -thickness_slope_kernel(n, thi)/2
#         j += 1
#     ind = indb[c]
#     b[i, 0] = airfoil.dydx[ind]
#     i += 1

# Trailing Edge Camber Fit
# j = 0
# for k in range(numc):
#     n = k+1
#     a[i, j] = camber_kernel(n, pi)
#     j += 1

# i = 0
# for c, xni in enumerate(xnt):
#     j = 0
#     for k in range(numc):
#         n = k+1
#         a[i, j] = camber_kernel(n, xni)
#         j += 1
#     for k in range(numt):
#         n = k
#         a[i, j] = thickness_kernel(n, xni)/2
#         j += 1
#     ind = indt[c]
#     b[i, 0] = airfoil.yn[ind]
#     i += 1
# for c, xni in enumerate(xnb):
#     j = 0
#     for k in range(numc):
#         n = k+1
#         a[i, j] = camber_kernel(n, xni)
#         j += 1
#     for k in range(numt):
#         n = k
#         a[i, j] = -thickness_kernel(n, xni)/2
#         j += 1
#     ind = indb[c]
#     b[i, 0] = airfoil.yn[ind]
#     i += 1

d = lstsq(a, b)[0]
# d = solve(a, b)
Cn = d.transpose().tolist()[0]

Bn = Cn[0:numc]

print(f'Bn = {Bn}')

An = Cn[numc:]

print(f'An = {An}')

#%% Plot Camber Line
nump = 500
thp = [pi*i/nump for i in range(nump+1)]
xnp = [0.5*(1-cos(thi)) for thi in thp]

ypc = []
mpc = []
ypt = []
ypb = []

for i in range(nump+1):
    thi = thp[i]
    xni = xnp[i]
    ypc.append(0.0)
    mpc.append(0.0)
    for j, Bnj in enumerate(Bn):
        n = j+1
        ypc[i] += Bnj*camber_kernel(n, thi)
        mpc[i] += Bnj*camber_slope_kernel(n, thi)
    ypt.append(ypc[i])
    ypb.append(ypc[i])    
    for j, Anj in enumerate(An):
        n = j+1
        ypt[i] += Anj*thickness_kernel(n, thi)/2
        ypb[i] -= Anj*thickness_kernel(n, thi)/2

print(f'mpc = {mpc}')

#%% Plot Profile
ax = airfoil.plot_normalised_aifoil()
ax.plot(xnp, ypc, label='Camber Line')
ax.plot(xnp, ypt, label='Top Surface')
ax.plot(xnp, ypb, label='Bottom Surface')
ax.scatter(xu, yu, label='Top Points')
ax.scatter(xl, yl, label='Bottom Points')
ax.set_ylim(-0.2, 0.2)
leg = ax.legend()
