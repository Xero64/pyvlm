from math import pi, sin, cos, acos, degrees
from numpy.linalg import lstsq
from numpy.matlib import zeros
from matplotlib.pyplot import figure
from pygeom.geom2d import CubicSpline2D, Point2D
from .spacing import full_cosine_spacing
from pygeom.geom1d import CubicSpline

class Airfoil(object):
    name = None
    x = None
    y = None
    xte = None
    yte = None
    xle = None
    yle = None
    ile = None
    crd = None
    xn = None
    yn = None
    spline = None
    splinec = None
    mle = None
    mte = None
    def __init__(self, name: str, x: list, y: list):
        self.name = name
        self.x = x
        self.y = y
        self.update()
    def update(self, num: int=80):
        self.xte = (self.x[0]+self.x[-1])/2
        self.yte = (self.y[0]+self.y[-1])/2
        self.xle = min(self.x)
        for i, yi in enumerate(self.y):
            if self.x[i] == self.xle:
                self.yle = yi
                self.ile = i
                break
        self.crd = self.xte-self.xle
        self.xn = [(xi-self.xle)/self.crd for xi in self.x]
        self.yn = [(yi-self.yle)/self.crd for yi in self.y]
        pnts = [Point2D(xi, yi) for xi, yi in zip(self.xn, self.yn)]
        self.spline = CubicSpline2D(pnts)
        dxle = self.spline.dr[self.ile].x
        dyle = self.spline.dr[self.ile].y
        self.mle = -dxle/dyle
        dx1 = self.spline.dr[0].x
        dy1 = self.spline.dr[0].y
        dx2 = self.spline.dr[-1].x
        dy2 = self.spline.dr[-1].y
        m1 = dy1/dx1
        m2 = dy2/dx2
        self.mte = (m1+m2)/2
        cspc = full_cosine_spacing(num)
        s = self.spline.arc_length()
        sle = s[self.ile]
        ste = s[-1]
        s1 = sle
        s2 = ste-sle
        slst1 = [cspci*s1 for cspci in cspc]
        slst2 = [cspci*s2+s1 for cspci in cspc]
        x1, y1 = self.spline.interpolate_spline_points(slst1)
        x2, y2 = self.spline.interpolate_spline_points(slst2)
        x1.reverse()
        y1.reverse()
        xc = [(x1i+x2i)/2 for x1i, x2i in zip(x1, x2)]
        yc = [(y1i+y2i)/2 for y1i, y2i in zip(y1, y2)]
        xnc, ync = [], []
        for xi, yi in zip(xc, yc):
            if 0.0 <= xi <= 1.0:
                xnc.append(xi)
                ync.append(yi)
        xnc[0], ync[0] = 0.0, 0.0
        self.splinec = CubicSpline(xnc, ync)
    def plot_airfoil(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
            ax.set_aspect('equal')
        ax.plot(self.x, self.y, label=self.name)
        return ax
    def plot_normalised_aifoil(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
            ax.set_aspect('equal')
        ax.plot(self.xn, self.yn, label=self.name)
        return ax
    def return_camber(self, xc: float):
        return self.splinec.single_interpolate_spline(xc)
    def return_camber_slope(self, xc: float):
        return self.splinec.single_interpolate_gradient(xc)
    def return_camber_angle(self, xc: float):
        from math import degrees, atan
        dydx = self.return_camber_slope(xc)
        return degrees(atan(dydx))

def airfoil_from_dat(datfilepath: str):
    x = []
    y = []
    with open(datfilepath, 'rt') as file:
        for i, line in enumerate(file):
            line = line.rstrip('\n')
            if i == 0:
                name = line.strip()
            else:
                split = line.split()
                if len(split) == 2:
                    x.append(float(split[0]))
                    y.append(float(split[1]))
    return Airfoil(name, x, y)
