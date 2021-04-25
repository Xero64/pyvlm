from pygeom.geom3d import Vector, ihat
from .latticecontrol import LatticeControl, latticecontrol_from_json
from ..tools.camber import FlatPlate, NACA4, NACA6Series
from ..tools import equal_spacing, full_cosine_spacing, semi_cosine_spacing
from ..tools.airfoil import airfoil_from_dat

class LatticeSection(object):
    pnt = None
    chord = None
    twist = None
    camber = None
    airfoil = None
    bspc = None
    yspace = None
    mirror = None
    noload = None
    ruled = None
    ctrls = None
    cdo = None
    bpos = None
    def __init__(self, pnt: Vector, chord: float, twist: float):
        self.pnt = pnt
        self.chord = chord
        self.twist = twist
        self.update()
    def update(self):
        self.noload = False
        self.mirror = False
        self.camber = FlatPlate()
        self.ctrls = {}
        self.cdo = 0.0
    def offset_position(self, xpos: float, ypos: float, zpos: float):
        self.pnt.x = self.pnt.x+xpos
        self.pnt.y = self.pnt.y+ypos
        self.pnt.z = self.pnt.z+zpos
    def offset_twist(self, twist: float):
        self.twist = self.twist+twist
    def set_span_equal_spacing(self, bnum: int):
        bsp = equal_spacing(2*bnum)
        self.bspc = [tuple(bsp[i*2:i*2+3]) for i in range(bnum)]
    def set_span_cosine_spacing(self, bnum: int):
        bsp = full_cosine_spacing(2*bnum)
        self.bspc = [tuple(bsp[i*2:i*2+3]) for i in range(bnum)]
    def set_span_semi_cosine_spacing(self, bnum: int):
        bsp = semi_cosine_spacing(2*bnum)
        self.bspc = [tuple(bsp[i*2:i*2+3]) for i in range(bnum)]
    def set_airfoil(self, airfoil: str):
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
    def set_noload(self, noload: bool):
        self.noload = noload
    def set_cdo(self, cdo: float):
        self.cdo = cdo
    def add_control(self, ctrl: LatticeControl):
        self.ctrls[ctrl.name] = ctrl
    def return_mirror(self):
        pnt = Vector(self.pnt.x, -self.pnt.y, self.pnt.z)
        chord = self.chord
        twist = self.twist
        sct = LatticeSection(pnt, chord, twist)
        sct.camber = self.camber
        sct.bspc = self.bspc
        sct.yspace = self.yspace
        sct.ctrls = self.ctrls
        sct.cdo = self.cdo
        sct.mirror = True
        return sct
    def return_point(self, percrd: float):
        return self.pnt+self.chord*percrd*ihat
    def get_camber(self, xc: float):
        return self.camber.cubic_interp(xc)
    def __repr__(self):
        return '<LatticeSection>'

def latticesection_from_json(sectdata: dict):
    if 'xpos' in sectdata:
        xpos = sectdata['xpos']
    else:
        return ValueError
    if 'ypos' in sectdata:
        ypos = sectdata['ypos']
    else:
        return ValueError
    if 'zpos' in sectdata:
        zpos = sectdata['zpos']
    else:
        return ValueError
    point = Vector(xpos, ypos, zpos)
    if 'chord' in sectdata:
        chord = sectdata['chord']
    else:
        return ValueError
    if 'twist' in sectdata:
        twist = sectdata['twist']
    else:
        twist = 0.0
    sct = LatticeSection(point, chord, twist)
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
