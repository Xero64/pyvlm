from math import asin, cos, pi
from pygeom.geom3d import Point, ihat
from pyvlm.tools.function import Function

class LatticeSection(object):
    pnt = None
    chord = None
    angle = None
    camber = None
    bspace = None
    yspace = None
    mirror = None
    def __init__(self, pnt: Point, chord: float, angle: float):
        self.pnt = pnt
        self.chord = chord
        self.angle = angle
        self.update()
    def update(self):
        self.mirror = False
        self.camber = Function([0.0, 1.0], [0.0, 0.0])
    def set_span_equal_spacing(self, numb: int):
        from pyvlm.tools import equal_spacing
        spc = equal_spacing(2*numb)
        self.bspace = spc[0::2]
        self.yspace = spc[1::2]
    def set_span_cosine_spacing(self, numb: int):
        from pyvlm.tools import full_cosine_spacing
        spc = full_cosine_spacing(2*numb)
        self.bspace = spc[0::2]
        self.yspace = spc[1::2]
    def set_span_semi_cosine_spacing(self, numb: int):
        from pyvlm.tools import semi_cosine_spacing
        spc = semi_cosine_spacing(2*numb)
        self.bspace = spc[0::2]
        self.yspace = spc[1::2]
    # def set_span_elliptical_spacing(self, numb: int):
    #     from pyvlm.tools import full_elliptical_spacing
    #     spc = full_elliptical_spacing(2*numb)
    #     self.bspace = spc[0::2]
    #     self.yspace = spc[1::2]
    # def set_span_semi_elliptical_spacing(self, numb: int):
    #     from pyvlm.tools import semi_elliptical_spacing
    #     spc = semi_elliptical_spacing(2*numb)
    #     self.bspace = spc[0::2]
    #     self.yspace = spc[1::2]
    def set_camber(self, xc: list, zc: list):
        xmin = min(xc)
        xmax = max(xc)
        xrng = xmax-xmin
        xnrm = [(xci-xmin)/xrng for xci in xc]
        znrm = [zci/xrng for zci in zc]
        self.camber = Function(xnrm, znrm)
    def return_mirror(self):
        pnt = Point(self.pnt.x, -self.pnt.y, self.pnt.z)
        chord = self.chord
        angle = self.angle
        sect = LatticeSection(pnt, chord, angle)
        sect.camber = self.camber
        sect.bspace = self.bspace
        sect.yspace = self.yspace
        sect.mirror = True
        return sect
    def return_point(self, percrd: float):
        return self.pnt+self.chord*percrd*ihat
    def get_camber(self, xc: float):
        return self.camber.cubic_interp(xc)
    def __repr__(self):
        return '<LatticeSection>'

def latticesecttion_from_json(sectdata: dict):
    xle = sectdata['xle']
    yle = sectdata['yle']
    zle = sectdata['zle']
    crd = sectdata['chord']
    if 'angle' in sectdata:
        ang = sectdata['angle']
    else:
        ang = 0.0
    pnt = Point(xle, yle, zle)
    sect = LatticeSection(pnt, crd, ang)
    if 'airfoil' in sectdata:
        sect.airfoil = sectdata['airfoil']
    if 'numb' in sectdata and 'bspace' in sectdata:
        numb = sectdata['numb']
        bspace = sectdata['bspace']
        if bspace == 'equal':
            sect.set_span_equal_spacing(numb)
        elif bspace == 'cosine':
            sect.set_span_cosine_spacing(numb)
        # elif bspace == 'elliptical':
        #     sect.set_span_elliptical_spacing(numb)
        elif bspace == 'semi-cosine':
            sect.set_span_semi_cosine_spacing(numb)
        # elif bspace == 'semi-elliptical':
        #     sect.set_span_semi_elliptical_spacing(numb)
    return sect
