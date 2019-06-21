from math import pi
from pygeom.geom3d import Point, Vector, ihat
from .latticepanel import LatticePanel

class LatticeStrip(object):
    lsid = None
    pnt1 = None
    pnt2 = None
    crd1 = None
    crd2 = None
    leni = None
    pnti = None
    pnls = None
    alpha = None
    msid = None
    sht = None
    def __init__(self, lsid: int, pnt1: Point, pnt2: Point, crd1: float, crd2: float):
        self.lsid = lsid
        self.pnt1 = pnt1
        self.pnt2 = pnt2
        self.crd1 = crd1
        self.crd2 = crd2
        self.update()
    def update(self):
        self.alpha = 0.0
        self.pnls = []
        self.leni = self.pnt2-self.pnt1
        self.pnti = self.pnt1+0.5*self.leni
        self.lent = Vector(0.0, self.leni.y, self.leni.z)
        self.dyt = self.leni.y
        self.dzt = self.leni.z
        self.dst = (self.dyt**2+self.dzt**2)**0.5
    def add_panel(self, pnl):
        pnl.strp = self
        self.pnls.append(pnl)
        self.sht.add_panel(pnl)
    def set_mirror(self, mstrp):
        self.msid = mstrp.lsid
        for i in range(len(self.pnls)):
            mpnl = mstrp.pnls[i]
            pnl = self.pnls[i]
            pnl.mpid = mpnl.lpid
    def trefftz_velocity(self, pnt: Point):
        r = Point(0.0, pnt.y, pnt.z)
        ra = Point(0.0, self.pnt1.y, self.pnt1.z)
        rb = Point(0.0, self.pnt2.y, self.pnt2.z)
        a = r-ra
        b = r-rb
        am = a.return_magnitude()
        bm = b.return_magnitude()
        x = ihat
        axx = a**x
        bxx = b**x
        vel = (axx/am**2-bxx/bm**2)/2/pi
        return vel
    @property
    def nrmt(self):
        return self.sht.cord.dirz
    def __repr__(self):
        return '<LatticeStrip {:d}>'.format(self.lsid)
    def __str__(self):
        frmstr = 'LatticeStrip\t{:d}\t{:.5f}\t{:.5f}\t{:}'
        return frmstr.format(self.lsid, self.pnt1, self.pnt2, self.msid)
    def __format__(self, format_spec: str):
        frmstr = 'LatticeStrip\t{:d}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:}'
        return frmstr.format(self.lsid, self.pnt1, self.pnt2, self.msid)
