from math import pi, atan2, degrees
from pygeom.geom3d import Point, Vector, ihat
from .latticepanel import LatticePanel

class LatticeStrip(object):
    lsid = None
    pnt1 = None
    pnt2 = None
    crd1 = None
    crd2 = None
    ang1 = None
    ang2 = None
    bspc = None
    leni = None
    pnti = None
    pnls = None
    msid = None
    sht = None
    lent = None
    dyt = None
    dzt = None
    dst = None
    _crd = None
    _ang = None
    _area = None
    _avecrd = None
    def __init__(self, lsid: int, pnt1: Point, pnt2: Point, crd1: float, crd2: float, bspc: float=0.5):
        self.lsid = lsid
        self.pnt1 = pnt1
        self.pnt2 = pnt2
        self.crd1 = crd1
        self.crd2 = crd2
        self.bspc = bspc
        self.update()
    def update(self):
        self.pnls = []
        self.leni = self.pnt2-self.pnt1
        self.pnti = self.pnt1+self.bspc*self.leni
        self.lent = Vector(0.0, self.leni.y, self.leni.z)
        self.dyt = self.leni.y
        self.dzt = self.leni.z
        self.dst = (self.dyt**2+self.dzt**2)**0.5
    def set_angles(self, ang1: float, ang2: float):
        self.ang1 = ang1
        self.ang2 = ang2
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
    def return_points(self, percrd: float):
        pnt1 = self.pnt1+self.crd1*percrd*ihat
        pnt2 = self.pnt2+self.crd2*percrd*ihat
        return pnt1, pnt2
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
    @property
    def dihedral(self):
        return degrees(atan2(self.dzt, self.dyt))
    @property
    def chord(self):
        if self._crd is None:
            self._crd = self.crd1+(self.crd2-self.crd1)*self.bspc
        return self._crd
    @property
    def avechord(self):
        if self._avecrd is None:
            self._avecrd = (self.crd1+self.crd2)/2
        return self._avecrd
    @property
    def angle(self):
        if self._ang is None:
            self._ang = self.ang1+(self.ang2-self.ang1)*self.bspc
        return self._ang
    @property
    def area(self):
        if self._area is None:
            self._area = self.dst*(self.crd1+self.crd2)/2
        return self._area
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
