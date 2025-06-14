from typing import TYPE_CHECKING, Any

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
    pnls: list['LatticePanel'] = None

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

    @classmethod
    def from_dict(cls, name: str, controldata: dict[str, Any]) -> 'LatticeControl':
        xhinge = controldata['xhinge']
        posgain = controldata.get('posgain', 1.0)
        neggain = controldata.get('neggain', 1.0)
        ctrl = cls(name, posgain, neggain, xhinge)
        hvec = Vector(0.0, 0.0, 0.0)
        if 'hvec' in controldata:
            hvec = Vector.from_dict(controldata['hvec'])
        ctrl.set_hinge_vector(hvec)
        ctrl.reverse = controldata.get('reverse', False)
        return ctrl
