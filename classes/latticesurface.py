from math import pi, cos
from .latticestrip import LatticeStrip
from .latticesheet import LatticeSheet
from .latticepanel import LatticePanel

class LatticeSurface(object):
    name = None
    sects = None
    cspace = None
    xspace = None
    bdsit = None
    strps = None
    pnts = None
    pnls = None
    msht = None
    msect = None
    mirrored = None
    def __init__(self, name: str, sects: list, mirror: bool):
        self.name = name
        self.sects = sects
        self.mirror = mirror
        self.update()
    def update(self):
        if self.mirror and self.sects[0].pnt.y == 0.0:
            numsect = len(self.sects)
            newsects = []
            for i in range(numsect-1):
                sect = self.sects[numsect-1-i]
                msect = sect.return_mirror()
                newsects.append(msect)
            for sect in self.sects:
                newsects.append(sect)
            self.sects = newsects
        elif self.mirror and self.sects[0].pnt.y != 0.0:
            print(f'Warning: Cannot mirror {self.name}.')
            self.mirror = False
    def set_chord_distribution(self, cspace: list):
        from pyvlm.tools import normalise_spacing
        self.cspace = normalise_spacing(cspace)
    def set_chord_equal_distribution(self, numc: int):
        from pyvlm.tools import equal_spacing
        self.cspace = equal_spacing(numc)
    def set_chord_cosine_distribution(self, numc: int):
        from pyvlm.tools import full_cosine_spacing
        if numc > 1:
            spc = full_cosine_spacing(4*numc+2)
            self.cspace = spc[1::4]
            self.cspace[0] = spc[0]
            self.cspace[-1] = spc[-1]
            self.xspace = spc[2::4]
        else:
            self.cspace = [0.0, 1.0]
            self.xspace = [0.25]
    def set_chord_elliptical_distribution(self, numc: int):
        from pyvlm.tools import full_elliptical_spacing
        self.cspace = full_elliptical_spacing(numc)
    def mesh(self, lsid: int, lpid: int):
        from pygeom.geom3d import Point
        from numpy.matlib import empty
        nums = len(self.sects)
        self.shts = []
        for i in range(nums-1):
            a, b = i, i+1
            secta = self.sects[a]
            sectb = self.sects[b]
            self.shts.append(LatticeSheet(secta, sectb))
        self.strps = []
        for sht in self.shts:
            lsid = sht.mesh_strips(lsid)
            self.strps += sht.strps
        pnts = [strp.pnt1 for strp in self.strps]
        pnts.append(self.strps[-1].pnt2)
        crds = [strp.crd1 for strp in self.strps]
        crds.append(self.strps[-1].crd2)
        lenb = len(pnts)
        lenc = len(self.cspace)
        self.pnts = empty((lenb, lenc), dtype=Point)
        for i in range(lenb):
            minx = pnts[i].x
            y = pnts[i].y
            z = pnts[i].z
            c = crds[i]
            for j in range(lenc):
                cd = self.cspace[j]
                x = minx+cd*c
                self.pnts[i, j] = Point(x, y, z)
        self.pnls = empty((lenb-1, lenc-1), dtype=LatticePanel)
        for i, strp in enumerate(self.strps):
            for j in range(lenc-1):
                pnts = [
                    self.pnts[i, j],
                    self.pnts[i+1, j],
                    self.pnts[i, j+1],
                    self.pnts[i+1, j+1]
                ]
                cspc = (self.xspace[j]-self.cspace[j])/(self.cspace[j+1]-self.cspace[j])
                pnl = LatticePanel(lpid, pnts, cspc=cspc, bspc=strp.bspc)
                self.pnls[i, j] = pnl
                lpid += 1
                strp.add_panel(pnl)
        if self.mirror:
            numstrp = len(self.strps)
            hlfstrp = int(numstrp/2)
            for i in range(hlfstrp):
                strp = self.strps[numstrp-1-i]
                mstrp = self.strps[i]
                strp.set_mirror(mstrp)
        for sht in self.shts:
            sht.set_control_points_and_normals()
        return lsid, lpid
    def point_xyz(self):
        from numpy.matlib import zeros
        x = zeros(self.pnts.shape)
        y = zeros(self.pnts.shape)
        z = zeros(self.pnts.shape)
        for i in range(self.pnts.shape[0]):
            for j in range(self.pnts.shape[1]):
                x[i, j] = self.pnts[i, j].x
                y[i, j] = self.pnts[i, j].y
                z[i, j] = self.pnts[i, j].z
        return x, y, z
    def return_panels(self):
        pnls = []
        shp = self.pnls.shape
        for i in range(shp[0]):
            for j in range(shp[1]):
                pnls.append(self.pnls[i, j])
        return pnls
    def print_geom(self):
        print(f'Printing {self.name:s} geometry.')
        print('Sections')
        for i, sect in enumerate(self.sects):
            print(f'{i:d}\t{sect.pnt:.5f}\t{sect.chord:.5f}')
        print('Sheets')
        for i, sht in enumerate(self.shts):
            print(f'{i:d}')
            print(f'{sht.sect1.pnt:.5f}\t{sht.sect1.chord:.5f}')
            print(f'{sht.sect2.pnt:.5f}\t{sht.sect2.chord:.5f}')
    def __repr__(self):
        return '<LatticeSurface {:s}>'.format(self.name)

def latticesurface_from_json(surfdata: dict, display: bool=False):
    from .latticesection import latticesecttion_from_json
    name = surfdata['name']
    if 'mirror' in surfdata:
        mirror = surfdata['mirror']
    else:
        mirror = False
    if display: print(f'Loading Surface: {name:s}')
    sects = []
    for sectdata in surfdata['sections']:
        sect = latticesecttion_from_json(sectdata)
        sects.append(sect)
    surf = LatticeSurface(name, sects, mirror)
    if 'numc' in surfdata and 'cspace' in surfdata:
        numc = surfdata['numc']
        cspace = surfdata['cspace']
        if cspace == 'equal':
            surf.set_chord_equal_distribution(numc)
        elif cspace == 'cosine':
            surf.set_chord_cosine_distribution(numc)
        elif cspace == 'elliptical':
            surf.set_chord_elliptical_distribution(numc)
    return surf
