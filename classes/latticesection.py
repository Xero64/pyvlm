from math import cos, pi
from pygeom.geom3d import Point
from pyvlm.tools.function import Function

class LatticeSection(object):
    pnt = None
    chord = None
    camber = None
    bdist = None
    def __init__(self, pnt: Point, chord: float):
        self.pnt = pnt
        self.chord = chord
        self.update()
    def update(self):
        self.camber = Function([0.0, 1.0], [0.0, 0.0])
    def set_span_distribution(self, cdist: list):
        self.bdist = cdist
    def set_span_sine_distribution(self, numb: int):
        self.bdist = [cos(i*pi/2/numb) for i in range(numb+1)]
        self.bdist.reverse()
    def set_span_equal_distribution(self, numb: int):
        self.bdist = [float(i)/numb for i in range(numb+1)]
    def set_span_cosine_distribution(self, numb: int):
        self.bdist = [0.5*(1.0-cos(i*pi/numb)) for i in range(numb+1)]
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
        sect = LatticeSection(pnt, chord)
        sect.camber = self.camber
        return sect
    def get_camber(self, xc: float):
        return self.camber.cubic_interp(xc)
    def __repr__(self):
        return '<LatticeSection>'

def latticesecttion_from_json(sectdata: dict):
    xle = sectdata['xle']
    yle = sectdata['yle']
    zle = sectdata['zle']
    crd = sectdata['chord']
    pnt = Point(xle, yle, zle)
    sect = LatticeSection(pnt, crd)
    if 'airfoil' in sectdata:
        sect.airfoil = sectdata['airfoil']
    if 'numb' in sectdata and 'bdist' in sectdata:
        numb = sectdata['numb']
        bdist = sectdata['bdist']
        if bdist == 'sine':
            sect.set_span_sine_distribution(numb)
        elif bdist == 'equal':
            sect.set_span_equal_distribution(numb)
        elif bdist == 'cosine':
            sect.set_span_cosine_distribution(numb)    
    return sect
