from typing import TYPE_CHECKING

from numpy import arctan2, degrees, pi, sqrt
from pygeom.geom3d import IHAT, Vector

if TYPE_CHECKING:
    from .latticepanel import LatticePanel
    from .latticesheet import LatticeSheet

class LatticeStrip():
    lsid: int = None
    pnt1: Vector = None
    pnt2: Vector = None
    crd1: float = None
    crd2: float = None
    ang1: float = None
    ang2: float = None
    cdo1: float = None
    cdo2: float = None
    bspc: tuple[float, float, float] = None
    bfrc: float = None
    leni: Vector = None
    pnti: Vector = None
    pntq: Vector = None
    pnls: list['LatticePanel'] = None
    msid: int = None
    sht: 'LatticeSheet' = None
    lent: Vector = None
    dyt: float = None
    dzt: float = None
    dst: float = None
    bpos: float = None
    _crd: float = None
    _ang: float = None
    _cdo: float = None
    _area: float = None

    def __init__(self, lsid: int, pnt1: Vector, pnt2: Vector, crd1: float,
                 crd2: float, bspc: tuple[float, float, float]) -> None:
        self.lsid = lsid
        self.pnt1 = pnt1
        self.pnt2 = pnt2
        self.crd1 = crd1
        self.crd2 = crd2
        self.bspc = bspc
        self.update()

    def update(self) -> None:
        self.bfrc = (self.bspc[1] - self.bspc[0])/(self.bspc[2] - self.bspc[0])
        self.pnls = []
        self.leni = self.pnt2 - self.pnt1
        self.pnti = self.pnt1 + self.bfrc*self.leni
        self.lent = Vector(0.0, self.leni.y, self.leni.z)
        self.dyt = self.leni.y
        self.dzt = self.leni.z
        self.dst = sqrt(self.dyt**2 + self.dzt**2)
        pnta = self.pnt1 + 0.25*self.crd1*IHAT
        pntb = self.pnt2 + 0.25*self.crd2*IHAT
        vecab = pntb - pnta
        self.pntq = pnta + self.bfrc*vecab

    def set_twists(self, ang1: float, ang2: float) -> None:
        self.ang1 = ang1
        self.ang2 = ang2

    def set_twist(self, twist: float) -> None:
        self._ang = twist

    def set_cdo(self, cdo1: float, cdo2: float) -> None:
        self.cdo1 = cdo1
        self.cdo2 = cdo2

    def add_panel(self, pnl) -> None:
        self.pnls.append(pnl)

    def return_points(self, percrd: float) -> tuple[Vector, Vector]:
        pnt1 = self.pnt1 + self.crd1*percrd*IHAT
        pnt2 = self.pnt2 + self.crd2*percrd*IHAT
        return pnt1, pnt2

    def trefftz_velocity(self, pnt: Vector) -> Vector:
        r = Vector(0.0, pnt.y, pnt.z)
        ra = Vector(0.0, self.pnt1.y, self.pnt1.z)
        rb = Vector(0.0, self.pnt2.y, self.pnt2.z)
        a = r-ra
        b = r-rb
        axx = Vector(0.0, a.z, -a.y)
        bxx = Vector(0.0, b.z, -b.y)
        am2 = a.dot(a)
        bm2 = b.dot(b)
        vel = (axx/am2-bxx/bm2)/2/pi
        return vel

    def trefftz_lift(self) -> float:
        if self.noload:
            return 0.0
        else:
            return self.lent.y

    def trefftz_yfrc(self) -> float:
        if self.noload:
            return 0.0
        else:
            return -self.lent.z

    def trefftz_drag(self, vel: float) -> float:
        if self.noload:
            return 0.0
        else:
            return -self.dst*vel/2

    @property
    def nrmt(self) -> Vector:
        return self.sht.cord.dirz

    @property
    def dihedral(self) -> float:
        return degrees(arctan2(self.dzt, self.dyt))

    @property
    def chord(self) -> float:
        if self._crd is None:
            self._crd = self.crd1 + (self.crd2 - self.crd1)*self.bfrc
        return self._crd

    @property
    def twist(self) -> float:
        if self._ang is None:
            self._ang = self.ang1 + (self.ang2 - self.ang1)*self.bfrc
        return self._ang

    @property
    def area(self) -> float:
        if self._area is None:
            self._area = self.dst*self.chord
        return self._area

    @property
    def cdo(self) -> float:
        if self._cdo is None:
            self._cdo = self.cdo1 + (self.cdo2 - self.cdo1)*self.bfrc
        return self._cdo

    @property
    def cdoarea(self) -> float:
        if self.noload:
            return 0.0
        else:
            return self.cdo*self.area

    @property
    def noload(self) -> bool:
        return self.sht.noload

    def __repr__(self) -> str:
        return '<LatticeStrip {:d}>'.format(self.lsid)

    def __str__(self) -> str:
        frmstr = 'LatticeStrip\t{:d}\t{:.5f}\t{:.5f}\t{:}'
        return frmstr.format(self.lsid, self.pnt1, self.pnt2, self.msid)

    def __format__(self, format_spec: str) -> float:
        frmstr = 'LatticeStrip\t{:d}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:'
        frmstr += format_spec
        frmstr += '}\t{:}'
        return frmstr.format(self.lsid, self.pnt1, self.pnt2, self.msid)
