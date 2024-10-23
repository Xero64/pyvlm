from matplotlib.pyplot import figure
from numpy import asarray, logical_and
from pygeom.geom1d import CubicSpline1D
from pygeom.geom2d import CubicSpline2D, Vector2D
from pygeom.tools.spacing import full_cosine_spacing


class Airfoil():
    name = None
    x = None
    y = None
    xte = None
    yte = None
    xpos = None
    ypos = None
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
        self.xpos = min(self.x)
        for i, yi in enumerate(self.y):
            if self.x[i] == self.xpos:
                self.ypos = yi
                self.ile = i
                break
        self.crd = self.xte-self.xpos
        self.xn = asarray([(xi-self.xpos)/self.crd for xi in self.x])
        self.yn = asarray([(yi-self.ypos)/self.crd for yi in self.y])
        # pnts = [Vector2D(xi, yi) for xi, yi in zip(self.xn, self.yn)]
        pnts = Vector2D(self.xn, self.yn)
        self.spline = CubicSpline2D(pnts)
        dr = self.spline.evaluate_first_derivatives_at_t(self.spline.s)
        dxpos = dr[self.ile].x
        dypos = dr[self.ile].y
        self.mle = -dxpos/dypos
        dx1 = dr[0].x
        dy1 = dr[0].y
        dx2 = dr[-1].x
        dy2 = dr[-1].y
        m1 = dy1/dx1
        m2 = dy2/dx2
        self.mte = (m1+m2)/2
        cspc = full_cosine_spacing(num)
        s = self.spline.s
        sle = s[self.ile]
        ste = s[-1]
        s1 = sle
        s2 = ste-sle
        slst1 = [cspci*s1 for cspci in cspc]
        slst2 = [cspci*s2+s1 for cspci in cspc]
        x1, y1 = self.spline.evaluate_points_at_t(slst1).to_xy()
        x2, y2 = self.spline.evaluate_points_at_t(slst2).to_xy()
        x1 = x1[::-1]
        y1 = y1[::-1]
        xc = (x1 + x2)/2
        yc = (y1 + y2)/2
        check = logical_and(xc >= 0.0, xc <= 1.0)
        xnc = xc[check]
        ync = yc[check]
        xnc[0], ync[0] = 0.0, 0.0
        self.splinec = CubicSpline1D(xnc, ync)

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
        return self.splinec.evaluate_points_at_t(xc)

    def return_camber_slope(self, xc: float):
        return self.splinec.evaluate_first_derivatives_at_t(xc)

    def return_camber_angle(self, xc: float):
        from math import atan, degrees
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
