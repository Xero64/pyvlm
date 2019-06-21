from math import pi, cos
from numpy import empty
from .latticesection import LatticeSection
from .latticepanel import LatticePanel
from .latticestrip import LatticeStrip
from pygeom.geom3d import ihat, Point, Coordinate

class LatticeSheet(object):
    sect1 = None
    sect2 = None
    bdist = None
    levec = None
    cord = None
    nrml = None
    strps = None
    pnls = None
    mind = None
    def __init__(self, sect1: LatticeSection, sect2: LatticeSection):
        self.sect1 = sect1
        self.sect2 = sect2
        self.update()
    def update(self):
        if self.sect1.bdist is not None:
            self.bdist = self.sect1.bdist
        elif self.sect2.bdist is not None:
            bdist = [1.0-bd for bd in self.sect2.bdist]
            bdist.reverse()
            self.bdist = bdist
        else:
            self.bdist = [0.0, 1.0]
        self.levec = self.sect2.pnt-self.sect1.pnt
        vecz = (ihat**self.levec).to_unit()
        vecy = (vecz**ihat).to_unit()
        self.cord = Coordinate(self.sect1.pnt, ihat, vecy, vecz)
        self.strps = []
        self.pnls = []
    def mesh_strips(self, lsid: int):
        self.strps = []
        pnta = self.sect1.pnt
        vecr = self.levec
        crda = self.sect1.chord
        crdb = self.sect2.chord
        crdr = crdb-crda
        lenb = len(self.bdist)
        for i in range(lenb-1):
            bd1 = self.bdist[i]
            bd2 = self.bdist[i+1]
            pnt1 = pnta+bd1*vecr
            pnt2 = pnta+bd2*vecr
            crd1 = crda+bd1*crdr
            crd2 = crda+bd2*crdr
            strp = LatticeStrip(lsid, pnt1, pnt2, crd1, crd2)
            strp.sht = self
            self.strps.append(strp)
            lsid += 1
        return lsid
    def add_panel(self, pnl: LatticePanel):
        pnl.sht = self
        self.pnls.append(pnl)
    def set_control_points_and_normals(self):
        for pnl in self.pnls:
            pntc = pnl.return_panel_point()
            # pnl.set_coordinate(self.cord)
            pnl.set_control_point(pntc)
            pnl.set_camber_angle(0.0)
    def __repr__(self):
        return '<LatticeSheet>'
