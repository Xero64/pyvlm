from math import pi, radians, cos, sin, atan2
from pygeom.geom3d import Point, Vector, Coordinate, ihat

class LatticePanel(object):
    lpid = None
    pnts = None
    bspc = None
    pntc = None
    pnti = None
    pnta = None
    pntb = None
    leni = None
    # lent = None
    # nrmt = None
    # dst = None
    # cord = None
    alc = None
    mpid = None
    strp = None
    sht = None
    crd = None
    area = None
    def __init__(self, lpid: int, pnts: list, bspc: float=0.5):
        self.lpid = lpid
        self.pnts = pnts
        self.bspc = bspc
        self.update()
    def update(self):
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        veca = pnt3-pnt1
        vecb = pnt4-pnt2
        self.pnta = pnt1+0.25*veca
        self.pntb = pnt2+0.25*vecb
        self.leni = self.pntb-self.pnta
        self.pnti = self.pnta+self.bspc*self.leni
        self.th = atan2(self.leni.z, self.leni.y)
        crda = veca.return_magnitude()
        crdb = vecb.return_magnitude()
        self.crd = crda+(crdb-crda)*self.bspc
        self.area = self.leni.return_magnitude()*(veca.return_magnitude()+vecb.return_magnitude())/2
        # self.lent = Vector(0.0, self.leni.y, self.leni.z)
        # self.nrmt = Vector(0.0, -self.leni.z, self.leni.y).to_unit()
        # self.nrmt = ihat**self.
        # self.dst = self.lent.return_magnitude()
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
    # def set_coordinate(self, cord: Coordinate):
    #     self.cord = cord
    def set_control_point(self, pntc: Point):
        self.pntc = pntc
    def set_camber_angle(self, alpha: float):
        self.alpha = alpha
    # def set_normal_vector(self, nrml: Vector):
    #     self.nrml = nrml.to_unit()
    @property
    def nrml(self):
        # alpha = self.alpha+self.strp.alpha
        # alrad = alpha*pi/180
        # alrad = radians(alpha)
        # cosal = cos(alrad)
        # sinal = sin(alrad)
        # lnrm = Vector(sinal, 0.0, cosal)
        # gnrm = self.cord.vector_to_global(lnrm)
        # return gnrm
        als = radians(self.strp.alpha)
        alc = radians(self.alpha)
        th = self.th
        return Vector(-sin(alc - als), -sin(th)*cos(alc - als), cos(th)*cos(alc - als))
    @property
    def dnda(self):
        # alpha = self.alpha+self.strp.alpha
        # # alrad = alpha*pi/180
        # alrad = radians(alpha)
        # cosal = cos(alrad)
        # sinal = sin(alrad)
        # ltan = Vector(cosal, 0.0, -sinal)
        # gtan = self.cord.vector_to_global(ltan)
        # return gtan
        als = radians(self.strp.alpha)
        alc = radians(self.alpha)
        th = self.th
        return Vector(cos(alc - als), -sin(th)*sin(alc - als), sin(alc - als)*cos(th))
    # def estimate_control_point(self):
    #     pnt1 = self.pnts[0]
    #     pnt2 = self.pnts[1]
    #     pnt3 = self.pnts[2]
    #     pnt4 = self.pnts[3]
    #     veca = pnt3-pnt1
    #     vecb = pnt4-pnt2
    #     pnta = pnt1+0.75*veca
    #     pntb = pnt2+0.75*vecb
    #     vecc = pntb-pnta
    #     self.pntc = pnta+0.5*vecc
    # def estimate_normal_vector(self):
    #     pnt1 = self.pnts[0]
    #     pnt2 = self.pnts[1]
    #     pnt3 = self.pnts[2]
    #     pnt4 = self.pnts[3]
    #     veca = pnt4-pnt1
    #     vecb = pnt3-pnt2
    #     self.nrml = (veca**vecb).to_unit()
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
    # def trefftz_velocity(self, pnt: Point):
    #     r = Point(0.0, pnt.y, pnt.z)
    #     ra = Point(0.0, self.pnta.y, self.pnta.z)
    #     rb = Point(0.0, self.pntb.y, self.pntb.z)
    #     a = r-ra
    #     b = r-rb
    #     am = a.return_magnitude()
    #     bm = b.return_magnitude()
    #     x = ihat
    #     axx = a**x
    #     bxx = b**x
    #     vel = (axx/am**2-bxx/bm**2)/2/pi
    #     return vel
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
