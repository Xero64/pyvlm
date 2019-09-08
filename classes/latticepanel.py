from math import pi, radians, cos, sin, atan2
from pygeom.geom3d import Point, Vector, Coordinate, ihat

class LatticePanel(object):
    lpid = None
    pnts = None
    cspc = None
    bspc = None
    pntc = None
    pnti = None
    pnta = None
    pntb = None
    leni = None
    alc = None
    mpid = None
    strp = None
    sht = None
    crd = None
    area = None
    def __init__(self, lpid: int, pnts: list, cspc: float=0.25, bspc: float=0.5):
        self.lpid = lpid
        self.pnts = pnts
        self.cspc = cspc
        self.bspc = bspc
        self.update()
    def update(self):
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        veca = pnt3-pnt1
        vecb = pnt4-pnt2
        self.pnta = pnt1+self.cspc*veca
        self.pntb = pnt2+self.cspc*vecb
        self.leni = self.pntb-self.pnta
        self.pnti = self.pnta+self.bspc*self.leni
        self.th = atan2(self.leni.z, self.leni.y)
        crda = veca.return_magnitude()
        crdb = vecb.return_magnitude()
        self.crd = crda+(crdb-crda)*0.5
        self.area = self.leni.return_magnitude()*(veca.return_magnitude()+vecb.return_magnitude())/2
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
    def set_control_point(self, pntc: Point):
        self.pntc = pntc
    def set_camber_angle(self, alpha: float):
        self.alpha = alpha
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
    def velocity(self, pnt: Point):
        r = pnt
        ra = self.pnta
        rb = self.pntb
        a = r-ra
        b = r-rb
        am = a.return_magnitude()
        bm = b.return_magnitude()
        x = ihat
        vel = Vector(0.0, 0.0, 0.0)
        if pnt != self.pnti:
            axb = a**b
            if axb.return_magnitude() != 0.0:
                den = am*bm+a*b
                vel += axb/den*(1/am+1/bm)
        axx = a**x
        if axx.return_magnitude() != 0.0:
            den = am-a*x
            vel += axx/den/am
        bxx = b**x
        if bxx.return_magnitude() != 0.0:
            den = bm-b*x
            vel -= bxx/den/bm
        vel = vel/(4*pi)
        return vel
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
