from typing import TYPE_CHECKING, Any, Dict, List, Tuple, Union

from pygeom.geom3d import IHAT, Vector

from ..tools import equal_spacing, full_cosine_spacing, semi_cosine_spacing
from ..tools.airfoil import airfoil_from_dat
from ..tools.camber import NACA4, FlatPlate, NACA6Series
from .latticecontrol import LatticeControl, latticecontrol_from_json

if TYPE_CHECKING:
    from ..tools.airfoil import Airfoil
    from .latticecontrol import LatticeControl
    AirfoilType = Union[Airfoil, FlatPlate, NACA4, NACA6Series]

class LatticeSection():
    pnt: Vector = None
    chord: float = None
    twist: float = None
    camber: 'AirfoilType' = None
    airfoil: str = None
    bspc: List[Tuple[float, float, float]] = None
    # yspace = None
    mirror: bool = None
    noload: bool = None
    ruled: bool = None
    ctrls: Dict[str, 'LatticeControl'] = None
    cdo: float = None
    bpos: float = None

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

    def set_airfoil(self, airfoil: str) -> None:
        if airfoil[-4:].lower() == '.dat':
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
        sct.bspc = self.bspc
        # sct.yspace = self.yspace
        sct.ctrls = self.ctrls
        sct.cdo = self.cdo
        sct.mirror = True
        return sct

    def return_point(self, percrd: float) -> Vector:
        return self.pnt + self.chord*percrd*IHAT

    # def get_camber(self, xc: float):
        # return self.camber.cubic_interp(xc)

    def __repr__(self):
        return '<LatticeSection>'

def latticesection_from_json(sectdata: Dict[str, Any]) -> LatticeSection:
    if 'xpos' in sectdata:
        xpos = sectdata['xpos']
    else:
        raise ValueError()
    if 'ypos' in sectdata:
        ypos = sectdata['ypos']
    else:
        raise ValueError()
    if 'zpos' in sectdata:
        zpos = sectdata['zpos']
    else:
        raise ValueError()
    point = Vector(xpos, ypos, zpos)
    if 'chord' in sectdata:
        chord = sectdata['chord']
    else:
        raise ValueError()
    if 'twist' in sectdata:
        twist = sectdata['twist']
    else:
        twist = 0.0
    sct = LatticeSection(point, chord, twist)
    sct.xoc = sectdata.get('xoc', 0.0)
    sct.zoc = sectdata.get('zoc', 0.0)


    if 'cdo' in sectdata:
        sct.set_cdo(sectdata['cdo'])
    if 'noload' in sectdata:
        sct.set_noload(sectdata['noload'])
    if 'airfoil' in sectdata:
        sct.set_airfoil(sectdata['airfoil'])
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
