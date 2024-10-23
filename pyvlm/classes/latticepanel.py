from typing import TYPE_CHECKING

from numpy import arctan2, cos, pi, radians, sin
from pygeom.geom3d import Vector

from .latticestrip import LatticeStrip

if TYPE_CHECKING:
    from .latticegrid import LatticeGrid
    from .latticesheet import LatticeSheet

TOL = 1e-12
FOURPI = 4*pi

class LatticePanel():
    lpid: int = None
    pnts: list['LatticeGrid'] = None
    cspc: list[float] = None
    cfr1: float = None
    cfr2: float = None
    pntc: Vector = None
    pnti: Vector = None
    pntg: Vector = None
    pnta: Vector = None
    pntb: Vector = None
    leni: Vector = None
    strp: LatticeStrip = None
    crd: float = None
    area: float = None

    def __init__(self, lpid: int, pnts: list['LatticeGrid'],
                 cspc: list[float], strp: LatticeStrip) -> None:
        self.lpid = lpid
        self.pnts = pnts
        self.cspc = cspc
        self.strp = strp
        self.update()

    def update(self) -> None:
        self.cfr1 = (self.cspc[1]-self.cspc[0])/(self.cspc[4]-self.cspc[0])
        self.cfr2 = (self.cspc[3]-self.cspc[0])/(self.cspc[4]-self.cspc[0])
        self.strp.add_panel(self)
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        veca = pnt3 - pnt1
        vecb = pnt4 - pnt2
        self.pnta = pnt1 + self.cfr1*veca
        self.pntb = pnt2 + self.cfr1*vecb
        self.leni = self.pntb - self.pnta
        self.pntg = self.pnta + 0.5*self.leni
        self.pnti = self.pnta + self.strp.bfrc*self.leni
        self.th = arctan2(self.leni.z, self.leni.y)
        self.width = (self.leni.y**2 + self.leni.z**2)**0.5
        crda = veca.return_magnitude()
        crdb = vecb.return_magnitude()
        self.crd = crda + (crdb - crda)*self.strp.bfrc
        self.area = self.strp.dst*self.crd
        al1 = self.sht.sct1.camber.return_camber_angle(self.cspc[3])
        al2 = self.sht.sct2.camber.return_camber_angle(self.cspc[3])
        self.alpha = al1 + self.strp.bspc[1]*(al2 - al1)
        pnta = pnt1 + self.cfr2*veca
        pntb = pnt2 + self.cfr2*vecb
        lenc = pntb - pnta
        self.pntc = pnta + self.strp.bfrc*lenc

    def return_panel_point(self) -> Vector:
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        veca = pnt3 - pnt1
        vecb = pnt4 - pnt2
        pnta = pnt1 + 0.75*veca
        pntb = pnt2 + 0.75*vecb
        lenc = pntb - pnta
        pntc = pnta + self.bspc*lenc
        return pntc

    @property
    def sht(self) -> 'LatticeSheet':
        return self.strp.sht

    @property
    def bspc(self) -> tuple[float, float, float]:
        return self.strp.bspc

    @property
    def nrml(self) -> Vector:
        als = radians(self.strp.twist)
        alc = radians(self.alpha)
        th = self.th
        return Vector(-sin(alc - als), -sin(th)*cos(alc - als), cos(th)*cos(alc - als))

    @property
    def dnda(self) -> Vector:
        als = radians(self.strp.twist)
        alc = radians(self.alpha)
        th = self.th
        return Vector(cos(alc - als), -sin(th)*sin(alc - als), sin(alc - als)*cos(th))

    @property
    def noload(self) -> bool:
        return self.strp.noload

    @property
    def cdoarea(self) -> float:
        if self.noload:
            return 0.0
        else:
            return self.strp.cdo*self.area

    def dndl(self, gain: float, hvec: Vector) -> Vector:
        return gain*hvec.cross(self.nrml)

    def velocity(self, pnt: Vector) -> Vector:
        vel = Vector(0.0, 0.0, 0.0)
        a = pnt - self.pnta
        b = pnt - self.pntb
        adb = a.dot(b)
        am = a.return_magnitude()
        bm = b.return_magnitude()
        denab = am*bm + adb
        if abs(denab) > TOL and am > TOL and bm > TOL:
            axb = a.cross(b)
            vel += axb/denab*(1/am + 1/bm)
        dena = am - a.x
        if abs(dena) > TOL:
            axx = Vector(0.0, a.z, -a.y)
            vel += axx/dena/am
        denb = bm - b.x
        if abs(denb) > TOL:
            bxx = Vector(0.0, b.z, -b.y)
            vel -= bxx/denb/bm
        vel = vel/FOURPI
        return vel

    def __repr__(self) -> str:
        return '<LatticePanel {:d}>'.format(self.lpid)

    def __str__(self) -> str:
        frmstr = 'LatticePanel\t{:d}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:}'
        pnt1 = self.pnts[0]
        pnt2 = self.pnts[1]
        pnt3 = self.pnts[2]
        pnt4 = self.pnts[3]
        return frmstr.format(self.lpid, pnt1, pnt2, pnt3, pnt4)

    def __format__(self, format_spec: str) -> str:
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
        return frmstr.format(self.lpid, pnt1, pnt2, pnt3, pnt4)
