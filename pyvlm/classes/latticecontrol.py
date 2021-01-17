from pygeom.geom3d import Vector

class LatticeControl(object):
    name = None
    posgain = None
    neggain = None
    xhinge = None
    uhvec = None
    reverse = None
    pnls = None
    def __init__(self, name: str, posgain: float, neggain: float, xhinge: float):
        self.name = name
        self.posgain = posgain
        self.neggain = neggain
        self.xhinge = xhinge
        self.pnls = []
    def set_hinge_vector(self, hvec: Vector):
        if hvec.return_magnitude() != 0.0:
            self.uhvec = hvec.to_unit()
        else:
            self.uhvec = hvec
    def duplicate(self, mirror=True):
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
    def add_panel(self, pnl):
        self.pnls.append(pnl)

def latticecontrol_from_json(name: str, controldata: dict):
    xhinge = controldata['xhinge']
    posgain = 1.0
    if 'posgain' in controldata:
        posgain = controldata['posgain']
    neggain = 1.0
    if 'neggain' in controldata:
        neggain = controldata['neggain']
    ctrl = LatticeControl(name, posgain, neggain, xhinge)
    hvec = Vector(0.0, 0.0, 0.0)
    if 'hvec' in controldata:
        hvec = vector_from_json(controldata['hvec'])
    ctrl.set_hinge_vector(hvec)
    reverse = False
    if 'reverse' in controldata:
        reverse = controldata['reverse']
    ctrl.reverse = reverse
    return ctrl

def vector_from_json(vectordata):
    x, y, z = None, None, None
    if 'x' in vectordata:
        x = vectordata['x']
    if 'y' in vectordata:
        y = vectordata['y']
    if 'z' in vectordata:
        z = vectordata['z']
    return Vector(x, y, z)
