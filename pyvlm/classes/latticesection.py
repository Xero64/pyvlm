from typing import TYPE_CHECKING, Any

from pygeom.geom3d import IHAT, Vector
from pygeom.tools.spacing import (equal_spacing, full_cosine_spacing,
                                  semi_cosine_spacing)

from ..tools.airfoil import airfoil_from_dat
from ..tools.camber import NACA4, FlatPlate, NACA6Series
from .latticecontrol import LatticeControl, latticecontrol_from_json

if TYPE_CHECKING:
    from ..tools.airfoil import Airfoil
    from .latticecontrol import LatticeControl

class LatticeSection():
    pnt: Vector = None
    chord: float = None
    twist: float = None
    camber: 'Airfoil | FlatPlate | NACA4 | NACA6Series' = None
    airfoil: str = None
    bspc: list[tuple[float, float, float]] = None
    mirror: bool = None
    noload: bool = None
    ruled: bool = None
    ctrls: dict[str, 'LatticeControl'] = None
    cdo: float = None
    bpos: float = None
    xoc: float = None
    zoc: float = None

    def __init__(self, pnt: Vector, chord: float, twist: float) -> None:
        self.pnt = pnt
        self.chord = chord
        self.twist = twist
        self.update()

    def update(self) -> None:
        self.noload = False
        self.mirror = False
        self.camber = FlatPlate()
        self.ctrls = {}
        self.cdo = 0.0

    def offset_position(self, xpos: float, ypos: float, zpos: float) -> None:
        self.pnt.x = self.pnt.x + xpos
        self.pnt.y = self.pnt.y + ypos
        self.pnt.z = self.pnt.z + zpos

    def offset_twist(self, twist: float) -> None:
        self.twist = self.twist + twist

    def set_span_equal_spacing(self, bnum: int) -> None:
        bsp = equal_spacing(2*bnum)
        self.bspc = [tuple(bsp[i*2:i*2+3]) for i in range(bnum)]

    def set_span_cosine_spacing(self, bnum: int) -> None:
        bsp = full_cosine_spacing(2*bnum)
        self.bspc = [tuple(bsp[i*2:i*2+3]) for i in range(bnum)]

    def set_span_semi_cosine_spacing(self, bnum: int) -> None:
        bsp = semi_cosine_spacing(2*bnum)
        self.bspc = [tuple(bsp[i*2:i*2+3]) for i in range(bnum)]

    def set_airfoil(self, airfoil: str | None) -> None:
        if airfoil is None:
            self.camber = FlatPlate()
            self.airfoil = None
        elif airfoil[-4:].lower() == '.dat':
            self.camber = airfoil_from_dat(airfoil)
            self.airfoil = airfoil
        elif airfoil[0:4].lower() == 'naca':
            code = airfoil[4:].strip()
            if len(code) == 4:
                self.camber = NACA4(code)
                self.airfoil = airfoil
            elif code[0] == '6':
                self.camber = NACA6Series(code)
                self.airfoil = airfoil

    def set_noload(self, noload: bool) -> None:
        self.noload = noload

    def set_cdo(self, cdo: float) -> None:
        self.cdo = cdo

    def add_control(self, ctrl: 'LatticeControl') -> None:
        self.ctrls[ctrl.name] = ctrl

    def return_mirror(self) -> 'LatticeSection':
        pnt = Vector(self.pnt.x, -self.pnt.y, self.pnt.z)
        chord = self.chord
        twist = self.twist
        sct = LatticeSection(pnt, chord, twist)
        sct.camber = self.camber
        sct.airfoil = self.airfoil
        sct.bspc = self.bspc
        sct.ctrls = self.ctrls
        sct.cdo = self.cdo
        sct.mirror = True
        sct.xoc = self.xoc
        sct.zoc = self.zoc
        return sct

    def return_point(self, percrd: float) -> Vector:
        return self.pnt + self.chord*percrd*IHAT

    # def get_camber(self, xc: float):
        # return self.camber.cubic_interp(xc)

    def __repr__(self):
        return '<LatticeSection>'

def latticesection_from_dict(sectdata: dict[str, Any],
                             defaults: dict[str, Any]) -> LatticeSection:
    """Create a LatticeSection object from a dictionary."""
    xpos = sectdata.get('xpos', None)
    ypos = sectdata.get('ypos', None)
    zpos = sectdata.get('zpos', None)
    point = Vector(xpos, ypos, zpos)
    chord = sectdata.get('chord', defaults.get('chord', None))
    twist = sectdata.get('twist', defaults.get('twist', None))
    sct = LatticeSection(point, chord, twist)
    sct.bpos = sectdata.get('bpos', None)
    sct.xoc = sectdata.get('xoc', defaults.get('xoc', None))
    sct.zoc = sectdata.get('zoc', defaults.get('zoc', None))
    sct.set_cdo(sectdata.get('cdo', defaults.get('cdo', 0.0)))
    sct.set_noload(sectdata.get('noload', defaults.get('noload', False)))
    sct.set_airfoil(sectdata.get('airfoil', defaults.get('airfoil', None)))
    if 'bnum' in sectdata and 'bspc' in sectdata:
        bnum = sectdata['bnum']
        bspc = sectdata['bspc']
        if bspc == 'equal':
            sct.set_span_equal_spacing(bnum)
        elif bspc in ('full-cosine', 'cosine'):
            sct.set_span_cosine_spacing(bnum)
        elif bspc == 'semi-cosine':
            sct.set_span_semi_cosine_spacing(bnum)
    if 'controls' in sectdata:
        for name in sectdata['controls']:
            ctrl = latticecontrol_from_json(name, sectdata['controls'][name])
            sct.add_control(ctrl)
    return sct
