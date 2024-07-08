from typing import TYPE_CHECKING, Any, Dict, List

from pygeom.geom3d import Vector

if TYPE_CHECKING:
    from .latticepanel import LatticePanel


class LatticeControl():
    name: str = None
    posgain: float = None
    neggain: float = None
    xhinge: float = None
    uhvec: Vector = None
    reverse: bool = None
    pnls: List['LatticePanel'] = None

    def __init__(self, name: str, posgain: float, neggain: float,
                 xhinge: float) -> None:
        self.name = name
        self.posgain = posgain
        self.neggain = neggain
        self.xhinge = xhinge
        self.pnls = []

    def set_hinge_vector(self, hvec: Vector) -> None:
        if hvec.return_magnitude() != 0.0:
            self.uhvec = hvec.to_unit()
        else:
            self.uhvec = hvec

    def duplicate(self, mirror: bool = True) -> 'LatticeControl':
        if mirror and self.reverse:
            posgain, neggain = -self.neggain, -self.posgain
        else:
            posgain, neggain = self.posgain, self.neggain
        if mirror:
            uhvec = Vector(self.uhvec.x, -self.uhvec.y, self.uhvec.z)
        else:
            uhvec = Vector(self.uhvec.x, self.uhvec.y, self.uhvec.z)
        ctrl = LatticeControl(self.name, posgain, neggain, self.xhinge)
        ctrl.reverse = False
        ctrl.set_hinge_vector(uhvec)
        return ctrl

    def add_panel(self, pnl: 'LatticePanel') -> None:
        self.pnls.append(pnl)

def latticecontrol_from_json(name: str, controldata: Dict[str, Any]) -> LatticeControl:
    xhinge = controldata['xhinge']
    posgain = controldata.get('posgain', 1.0)
    neggain = controldata.get('neggain', 1.0)
    ctrl = LatticeControl(name, posgain, neggain, xhinge)
    hvec = Vector(0.0, 0.0, 0.0)
    if 'hvec' in controldata:
        hvec = vector_from_json(controldata['hvec'])
    ctrl.set_hinge_vector(hvec)
    ctrl.reverse = controldata.get('reverse', False)
    return ctrl

def vector_from_json(vectordata: Dict[str, float]) -> Vector:
    x, y, z = None, None, None
    if 'x' in vectordata:
        x = vectordata['x']
    if 'y' in vectordata:
        y = vectordata['y']
    if 'z' in vectordata:
        z = vectordata['z']
    return Vector(x, y, z)
