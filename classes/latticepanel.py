from .latticestrip import LatticeStrip
from math import pi, radians, cos, sin, atan2
from pygeom.geom3d import Point, Vector, Coordinate, ihat

class LatticePanel(object):
    lpid = None
    pnts = None
    cspc = None
    cfr1 = None
    cfr2 = None
    pntc = None
    pnti = None
    pntg = None
    pnta = None
    pntb = None
    leni = None
    slope = None
    mpid = None
    strp = None
    crd = None
    area = None
    def __init__(self, lpid: int, pnts: list, cspc: float, strp: LatticeStrip):
        self.lpid = lpid
        self.pnts = pnts
        self.cspc = cspc
        self.strp = strp
        self.update()
    def update(self):
        self.cfra = (self.cspc[1]-self.cspc[0])/(self.cspc[4]-self.cspc[0])
        self.cfrb = (self.cspc[3]-self.cspc[0])/(self.cspc[4]-self.cspc[0])
        self.strp.add_panel(self)
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        veca = pnt3-pnt1
        vecb = pnt4-pnt2
        self.pnta = pnt1+self.cfra*veca
        self.pntb = pnt2+self.cfra*vecb
        self.leni = self.pntb-self.pnta
        self.pntg = self.pnta+0.5*self.leni
        self.pnti = self.pnta+self.strp.bfrc*self.leni
        self.th = atan2(self.leni.z, self.leni.y)
        self.width = (self.leni.y**2+self.leni.z**2)**0.5
        crda = veca.return_magnitude()
        crdb = vecb.return_magnitude()
        self.crd = crda+(crdb-crda)*self.strp.bfrc
        self.area = self.strp.dst*self.crd
        # self.area = self.leni.return_magnitude()*(veca.return_magnitude()+vecb.return_magnitude())/2
        al1 = self.sht.sect1.camber.return_camber_angle(self.cspc[3])
        al2 = self.sht.sect2.camber.return_camber_angle(self.cspc[3])
        self.alpha = self.strp.bspc[1]*(al2-al1)+al1
        pnta = pnt1+self.cfrb*veca
        pntb = pnt2+self.cfrb*vecb
        lenc = pntb-pnta
        self.pntc = pnta+self.strp.bfrc*lenc
    def return_panel_point(self):
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        veca = pnt3-pnt1
        vecb = pnt4-pnt2
        pnta = pnt1+0.75*veca
        pntb = pnt2+0.75*vecb
        lenc = pntb-pnta
        pntc = pnta+self.bspc*lenc
        return pntc
    @property
    def sht(self):
        return self.strp.sht
    @property
    def bspc(self):
        return self.strp.bspc
    @property
    def yspc(self):
        return self.strp.yspc
    @property
    def nrml(self):
        als = radians(self.strp.angle)
        alc = radians(self.alpha)
        th = self.th
        return Vector(-sin(alc - als), -sin(th)*cos(alc - als), cos(th)*cos(alc - als))
    @property
    def dnda(self):
        als = radians(self.strp.angle)
        alc = radians(self.alpha)
        th = self.th
        return Vector(cos(alc - als), -sin(th)*sin(alc - als), sin(alc - als)*cos(th))
    @property
    def noload(self):
        return self.strp.noload
    def dndl(self, gain: float, hvec: Vector):
        return gain*(hvec**self.nrml)
    def velocity(self, pnt: Point):
        r = pnt
        ra = self.pnta
        rb = self.pntb
        a = r-ra
        b = r-rb
        am = a.return_magnitude()
        bm = b.return_magnitude()
        # x = ihat
        vel = Vector(0.0, 0.0, 0.0)
        if pnt != self.pnti:
            axb = a**b
            if axb.return_magnitude() != 0.0:
                den = am*bm+a*b
                vel += axb/den*(1/am+1/bm)
        # adx = a*x
        # axx = a**x
        adx = a.x
        axx = Vector(0.0, a.z, -a.y)
        if axx.return_magnitude() != 0.0:
            den = am-adx
            vel += axx/den/am
        # bdx = b*x
        # bxx = b**x
        bdx = b.x
        bxx = Vector(0.0, b.z, -b.y)
        if bxx.return_magnitude() != 0.0:
            den = bm-bdx
            vel -= bxx/den/bm
        vel = vel/(4*pi)
        return vel
    def induced_force(self, vel: Vector):
        if self.noload:
            return Vector(0.0, 0.0, 0.0)
        else:
            return vel**self.leni
    def induced_velocity(self, pnt: Point):
        if self.noload:
            return Vector(0.0, 0.0, 0.0)
        else:
            return self.velocity(pnt)
    def __repr__(self):
        return '<LatticePanel {:d}>'.format(self.lpid)
    def __str__(self):
        frmstr = 'LatticePanel\t{:d}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:}'
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        return frmstr.format(self.lpid, pnt1, pnt2, pnt3, pnt4, self.mpid)
    def __format__(self, format_spec: str):
        frmstr = 'LatticePanel\t{:d}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:}'
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        return frmstr.format(self.lpid, pnt1, pnt2, pnt3, pnt4, self.mpid)
