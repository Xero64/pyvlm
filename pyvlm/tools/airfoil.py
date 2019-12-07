from math import atan2, tan, pi, cos
from matplotlib.pyplot import figure
from pygeom.geom2d import CubicSpline2D, Point2D, InfLine2D, Vector2D

class Airfoil(object):
    name = None
    x = None
    y = None
    xte = None
    yte = None
    xle = None
    yle = None
    ile = None
    # slp = None
    crd = None
    xn = None
    yn = None
    # rngt = None
    # rngb = None
    mte = None
    spline = None
    mte = None
    # dydx = None
    # xu = None
    # yu = None
    # xl = None
    # yl = None
    def __init__(self, name: str, x: list, y: list):
        self.name = name
        self.x = x
        self.y = y
        self.update()
    def update(self):
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
        # self.slp = (yn[0]+yn[-1])/2
        # self.yn = [yn[i]-xi*self.slp for i, xi in enumerate(self.xn)]

        # num = len(self.xn)
        # area1 = 0.0
        # for i in range(0, self.ile):
        #     area1 -= (self.yn[i+1]+self.yn[i])/2*(self.xn[i+1]-self.xn[i])
        # area2 = 0.0
        # for i in range(self.ile, num-1):
        #     area2 -= (self.yn[i+1]+self.yn[i])/2*(self.xn[i+1]-self.xn[i])
        # if area1 > area2:
        #     self.rngt = range(1, self.ile)
        #     self.rngb = range(self.ile+1, num-1)
        # else:
        #     self.rngb = range(1, self.ile)
        #     self.rngt = range(self.ile+1, num-1)

        pnts = [Point2D(self.xn[i], self.yn[i]) for i in range(len(self.xn))]
        self.spline = CubicSpline2D(pnts)

        dx1 = self.spline.dr[0].x
        dy1 = self.spline.dr[0].y

        dx2 = self.spline.dr[-1].x
        dy2 = self.spline.dr[-1].y
        
        # dx1 = self.xn[-1]-self.xn[-2]
        # dy1 = self.yn[-1]-self.yn[-2]
        # s1 = (dx1**2+dy1**2)**0.5
        # dx1 = dx1/s1
        # dy1 = dy1/s1
        # dx2 = self.xn[0]-self.xn[1]
        # dy2 = self.yn[0]-self.yn[1]
        # s2 = (dx2**2+dy2**2)**0.5
        # dx2 = dx2/s2
        # dy2 = dy2/s2
        self.mte = (dy1/dx1+dy2/dx2)/2
        
        # print(f'th = {th}')
    def cosine_spline(self, num: int=40):
        # self.dydx = [self.spline.dr[i].y/self.spline.dr[i].x for i in range(num)]
        th = [i*pi/num for i in range(num+1)]
        xu = []
        xl = []
        yu = []
        yl = []
        for i in range(1, num):
            thi = th[i]
            # print(f'\nthi = {thi}')
            xi = 0.5*(1-cos(thi))
            # print(f'xi = {xi}')
            infln = InfLine2D(Point2D(xi, 0.0), Vector2D(0.0, 1.0))
            pnts = self.spline.line_intersection_points(infln)
            # print(len(pnts))
            # print(pnts)
            if len(pnts) == 2:
                # if abs(pnts[0].y - pnts[1].y) < 1e-6:
                #     pass
                if pnts[0].y > pnts[1].y:
                    pntu = pnts[0]
                    pntl = pnts[1]
                else:
                    pntu = pnts[1]
                    pntl = pnts[0]
                xu.append(pntu.x)
                yu.append(pntu.y)
                xl.append(pntl.x)
                yl.append(pntl.y)
        return xu, yu, xl, yl
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
