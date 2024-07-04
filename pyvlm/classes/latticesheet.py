from typing import TYPE_CHECKING, List, Dict, Tuple
from numpy import arctan2, cos, degrees, radians, sin

from pygeom.geom3d import IHAT, Vector, Coordinate

from .latticestrip import LatticeStrip

if TYPE_CHECKING:
    from .latticesection import LatticeSection
    from .latticecontrol import LatticeControl
    from .latticepanel import LatticePanel


class LatticeSheet():
    sct1: 'LatticeSection' = None
    sct2: 'LatticeSection' = None
    bspc: List[Tuple[float, float, float]] = None
    mirror: bool = None
    levec: Vector = None
    cord: Coordinate = None
    strps: List['LatticeStrip'] = None
    pnls: List['LatticePanel'] = None
    ctrls: Dict[str, 'LatticeControl'] = None
    noload: bool = None
    ruled: bool = None
    width: float = None
    area: float = None

    def __init__(self, sct1: 'LatticeSection', sct2: 'LatticeSection') -> None:
        self.sct1 = sct1
        self.sct2 = sct2
        self.update()

    def update(self) -> None:
        if self.sct1.mirror or self.sct2.mirror:
            self.mirror = True
        else:
            self.mirror = False
        self.inherit_ruled()
        self.inherit_noload()
        self.inherit_spacing()
        self.inherit_controls()
        self.levec = self.sct2.pnt - self.sct1.pnt
        vecz = IHAT.cross(self.levec).to_unit()
        vecy = vecz.cross(IHAT).to_unit()
        self.cord = Coordinate(self.sct1.pnt, IHAT, vecy)
        self.strps = []
        self.pnls = []
        self.width = self.levec.dot(self.cord.diry)
        self.area = self.width*(self.sct2.chord+self.sct1.chord)/2

    def mesh_strips(self, lsid: int) -> None:
        self.strps = []
        pnta = self.sct1.pnt
        vecr = self.levec
        crda = self.sct1.chord
        crdb = self.sct2.chord
        crdr = crdb - crda
        anga = self.sct1.twist
        angb = self.sct2.twist
        angr = angb - anga
        if self.ruled:
            radanga = radians(anga)
            cosanga = cos(radanga)
            sinanga = sin(radanga)
            xla = cosanga*crda
            zla = -sinanga*crda
            radangb = radians(angb)
            cosangb = cos(radangb)
            sinangb = sin(radangb)
            xlb = cosangb*crdb
            zlb = -sinangb*crdb
            xlr = xlb - xla
            zlr = zlb - zla
        cdoa = self.sct1.cdo
        cdob = self.sct2.cdo
        xoca = self.sct1.xoc
        xocb = self.sct2.xoc
        zoca = self.sct1.zoc
        zocb = self.sct2.zoc
        cdor = cdob - cdoa
        xocr = xocb - xoca
        zocr = zocb - zoca
        lenb = len(self.bspc)
        for i in range(lenb):
            bspc = self.bspc[i]
            bsp1 = bspc[0]
            bspm = bspc[1]
            bsp2 = bspc[2]
            xoc1 = xoca + bsp1*xocr
            xoc2 = xoca + bsp2*xocr
            zoc1 = zoca + bsp1*zocr
            zoc2 = zoca + bsp2*zocr
            crd1 = crda + bsp1*crdr
            crd2 = crda + bsp2*crdr
            pnt1 = pnta + bsp1*vecr - Vector(xoc1, 0.0, zoc1)*crd1
            pnt2 = pnta + bsp2*vecr - Vector(xoc2, 0.0, zoc2)*crd2
            strp = LatticeStrip(lsid, pnt1, pnt2, crd1, crd2, bspc)
            if self.ruled:
                xl1 = xla + bsp1*xlr
                xl2 = xla + bsp2*xlr
                xlm = xla + bspm*xlr
                zl1 = zla + bsp1*zlr
                zl2 = zla + bsp2*zlr
                zlm = zla + bspm*zlr
                ang1 = degrees(arctan2(-zl1, xl1))
                ang2 = degrees(arctan2(-zl2, xl2))
                angm = degrees(arctan2(-zlm, xlm))
            else:
                ang1 = anga + bsp1*angr
                ang2 = anga + bsp2*angr
                angm = anga + bspm*angr
            strp.set_twists(ang1, ang2)
            strp.set_twist(angm)
            cdo1 = cdoa + bsp1*cdor
            cdo2 = cdoa + bsp2*cdor
            strp.set_cdo(cdo1, cdo2)
            strp.sht = self
            self.strps.append(strp)
            lsid += 1
        return lsid

    def inherit_panels(self) -> None:
        self.pnls = []
        for strp in self.strps:
            for pnl in strp.pnls:
                self.pnls.append(pnl)

    def inherit_ruled(self) -> None:
        if self.mirror:
            self.ruled = self.sct2.ruled
        else:
            self.ruled = self.sct1.ruled

    def inherit_noload(self) -> None:
        if self.mirror:
            self.noload = self.sct2.noload
        else:
            self.noload = self.sct1.noload

    def inherit_spacing(self) -> None:
        if self.mirror:
            if self.sct2.bspc is None:
                self.bspc = [(0.0, 0.5, 1.0)]
            else:
                lenb = len(self.sct2.bspc)
                self.bspc = []
                for i in range(lenb):
                    blst = []
                    for j in range(2, -1, -1):
                        blst.append(1.0 - self.sct2.bspc[i][j])
                    self.bspc.append(tuple(blst))
                self.bspc.reverse()
        else:
            if self.sct1.bspc is None:
                self.bspc = [(0.0, 0.5, 1.0)]
            else:
                self.bspc = self.sct1.bspc

    def inherit_controls(self) -> None:
        self.ctrls = {}
        if self.mirror:
            for control in self.sct2.ctrls:
                ctrl = self.sct2.ctrls[control]
                newctrl = ctrl.duplicate(mirror=True)
                self.ctrls[control] = newctrl
        else:
            for control in self.sct1.ctrls:
                ctrl = self.sct1.ctrls[control]
                newctrl = ctrl.duplicate(mirror=False)
                self.ctrls[control] = newctrl
        for control in self.ctrls:
            ctrl = self.ctrls[control]
            if ctrl.uhvec.return_magnitude() == 0.0:
                pnt1 = self.sct1.pnt
                crd1 = self.sct1.chord
                pnta = pnt1 + crd1*IHAT*ctrl.xhinge
                pnt2 = self.sct2.pnt
                crd2 = self.sct2.chord
                pntb = pnt2 + crd2*IHAT*ctrl.xhinge
                hvec = pntb - pnta
                ctrl.set_hinge_vector(hvec)

    def set_control_panels(self) -> None:
        for control in self.ctrls:
            ctrl = self.ctrls[control]
            for pnl in self.pnls:
                if pnl.cspc[3] >= ctrl.xhinge:
                    ctrl.add_panel(pnl)

    def set_strip_bpos(self) -> None:
        bpos = self.sct1.bpos
        for i, strp in enumerate(self.strps):
            bspc = self.bspc[i]
            strp.bpos = bpos + self.width*bspc[1]

    def __repr__(self):
        return '<LatticeSheet>'
