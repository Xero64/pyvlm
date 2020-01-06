from math import asin, cos, pi
from pygeom.geom3d import Point, ihat
from pyvlm.tools.camber import FlatPlate
from .latticecontrol import LatticeControl, latticecontrol_from_json

class LatticeSection(object):
    pnt = None
    chord = None
    angle = None
    camber = None
    bspace = None
    yspace = None
    mirror = None
    noload = None
    ctrls = None
    cdo = None
    bpos = None
    def __init__(self, pnt: Point, chord: float, angle: float):
        self.pnt = pnt
        self.chord = chord
        self.angle = angle
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
    def offset_angle(self, angle: float):
        self.angle = self.angle+angle
    def set_span_equal_spacing(self, numb: int):
        from pyvlm.tools import equal_spacing
        bsp = equal_spacing(2*numb)
        self.bspace = [tuple(bsp[i*2:i*2+3]) for i in range(numb)]
    def set_span_cosine_spacing(self, numb: int):
        from pyvlm.tools import full_cosine_spacing
        bsp = full_cosine_spacing(2*numb)
        self.bspace = [tuple(bsp[i*2:i*2+3]) for i in range(numb)]
    def set_span_semi_cosine_spacing(self, numb: int):
        from pyvlm.tools import semi_cosine_spacing
        bsp = semi_cosine_spacing(2*numb)
        self.bspace = [tuple(bsp[i*2:i*2+3]) for i in range(numb)]
    def set_airfoil(self, airfoil: str):
        if airfoil[0:4].lower() == 'naca':
            code = airfoil[4:].strip()
            if len(code) == 4:
                from ..tools.camber import NACA4
                self.camber = NACA4(code)
                self.airfoil = airfoil
            elif code[0] == '6':
                from ..tools.camber import NACA6Series
                self.camber = NACA6Series(code)
                self.airfoil = airfoil
        elif airfoil[-4:].lower() == '.dat':
            from ..tools.airfoil import airfoil_from_dat
            self.camber = airfoil_from_dat(airfoil)
            self.airfoil = airfoil
    def set_noload(self, noload: bool):
        self.noload = noload
    def set_cdo(self, cdo: float):
        self.cdo = cdo
    def add_control(self, ctrl: LatticeControl):
        self.ctrls[ctrl.name] = ctrl
    def return_mirror(self):
        pnt = Point(self.pnt.x, -self.pnt.y, self.pnt.z)
        chord = self.chord
        angle = self.angle
        sect = LatticeSection(pnt, chord, angle)
        sect.camber = self.camber
        sect.bspace = self.bspace
        sect.yspace = self.yspace
        sect.ctrls = self.ctrls
        sect.cdo = self.cdo
        sect.mirror = True
        return sect
    def return_point(self, percrd: float):
        return self.pnt+self.chord*percrd*ihat
    def get_camber(self, xc: float):
        return self.camber.cubic_interp(xc)
    def __repr__(self):
        return '<LatticeSection>'

def latticesecttion_from_json(sectdata: dict):
    if 'xle' in sectdata:
        xle = sectdata['xle']
    else:
        return ValueError
    if 'yle' in sectdata:
        yle = sectdata['yle']
    else:
        return ValueError
    if 'zle' in sectdata:
        zle = sectdata['zle']
    else:
        return ValueError
    point = Point(xle, yle, zle)
    if 'chord' in sectdata:
        chord = sectdata['chord']
    else:
        return ValueError
    if 'angle' in sectdata:
        angle = sectdata['angle']
    else:
        angle = 0.0
    sect = LatticeSection(point, chord, angle)
    if 'cdo' in sectdata:
        sect.set_cdo(sectdata['cdo'])
    if 'noload' in sectdata:
        sect.set_noload(sectdata['noload'])
    if 'airfoil' in sectdata:
        sect.set_airfoil(sectdata['airfoil'])
    if 'numb' in sectdata and 'bspace' in sectdata:
        numb = sectdata['numb']
        bspace = sectdata['bspace']
        if bspace == 'equal':
            sect.set_span_equal_spacing(numb)
        elif bspace == 'cosine':
            sect.set_span_cosine_spacing(numb)
        elif bspace == 'semi-cosine':
            sect.set_span_semi_cosine_spacing(numb)
    if 'controls' in sectdata:
        for name in sectdata['controls']:
            ctrl = latticecontrol_from_json(name, sectdata['controls'][name])
            sect.add_control(ctrl)
    return sect
