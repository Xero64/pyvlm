from math import pi, atan2, degrees
from pygeom.geom3d import Point, Vector, ihat

class LatticeStrip(object):
    lsid = None
    pnt1 = None
    pnt2 = None
    crd1 = None
    crd2 = None
    ang1 = None
    ang2 = None
    cdo1 = None
    cdo2 = None
    bspc = None
    bfrc = None
    leni = None
    pnti = None
    pntq = None
    # pnta = None
    # pntb = None
    # teni = None
    # pnte = None
    pnls = None
    msid = None
    sht = None
    lent = None
    dyt = None
    dzt = None
    dst = None
    bpos = None
    _crd = None
    _ang = None
    _cdo = None
    _area = None
    def __init__(self, lsid: int, pnt1: Point, pnt2: Point, crd1: float, crd2: float, bspc: tuple):
        self.lsid = lsid
        self.pnt1 = pnt1
        self.pnt2 = pnt2
        self.crd1 = crd1
        self.crd2 = crd2
        self.bspc = bspc
        self.update()
    def update(self):
        self.bfrc = (self.bspc[1]-self.bspc[0])/(self.bspc[2]-self.bspc[0])
        self.pnls = []
        self.leni = self.pnt2-self.pnt1
        self.pnti = self.pnt1+self.bfrc*self.leni
        self.lent = Vector(0.0, self.leni.y, self.leni.z)
        self.dyt = self.leni.y
        self.dzt = self.leni.z
        self.dst = (self.dyt**2+self.dzt**2)**0.5
        pnta = self.pnt1+0.25*self.crd1*ihat
        pntb = self.pnt2+0.25*self.crd2*ihat
        vecab = pntb-pnta
        self.pntq = pnta+self.bfrc*vecab
        # self.pnta = self.pnt1+self.crd1*ihat
        # self.pntb = self.pnt2+self.crd2*ihat
        # self.teni = self.pntb-self.pnta
        # self.pnte = self.pnta+self.bfrc*self.teni+1.0*self.chord*ihat
    def set_angles(self, ang1: float, ang2: float):
        self.ang1 = ang1
        self.ang2 = ang2
    def set_cdo(self, cdo1: float, cdo2: float):
        self.cdo1 = cdo1
        self.cdo2 = cdo2
    def add_panel(self, pnl):
        self.pnls.append(pnl)
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
        axx = Vector(0.0, a.z, -a.y)
        bxx = Vector(0.0, b.z, -b.y)
        am2 = a*a
        bm2 = b*b
        vel = (axx/am2-bxx/bm2)/2/pi
        return vel
    def trefftz_lift(self):
        if self.noload:
            return 0.0
        else:
            return self.lent.y
    def trefftz_yfrc(self):
        if self.noload:
            return 0.0
        else:
            return -self.lent.z
    def trefftz_drag(self, vel: float):
        if self.noload:
            return 0.0
        else:
            return -self.dst*vel/2
    @property
    def nrmt(self):
        return self.sht.cord.dirz
    @property
    def dihedral(self):
        return degrees(atan2(self.dzt, self.dyt))
    @property
    def chord(self):
        if self._crd is None:
            self._crd = self.crd1+(self.crd2-self.crd1)*self.bfrc
        return self._crd
    @property
    def angle(self):
        if self._ang is None:
            self._ang = self.ang1+(self.ang2-self.ang1)*self.bfrc
        return self._ang
    @property
    def area(self):
        if self._area is None:
            self._area = self.dst*self.chord
        return self._area
    @property
    def cdo(self):
        if self._cdo is None:
            self._cdo = self.cdo1+(self.cdo2-self.cdo1)*self.bfrc
        return self._cdo
    @property
    def cdoarea(self):
        if self.noload:
            return 0.0
        else:
            return self.cdo*self.area
    @property
    def noload(self):
        return self.sht.noload
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
