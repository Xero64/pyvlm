from math import pi, cos
from numpy import empty
from .latticesection import LatticeSection
from .latticepanel import LatticePanel
from .latticestrip import LatticeStrip
from pygeom.geom3d import ihat, Point, Coordinate

class LatticeSheet(object):
    sect1 = None
    sect2 = None
    bspace = None
    yspace = None
    mirror = None
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
        if self.sect1.mirror or self.sect2.mirror:
            self.mirror = True
        else:
            self.mirror = False
        self.bspace = inherit_spacing(self)
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
        anga = self.sect1.angle
        angb = self.sect2.angle
        angr = angb-anga
        lenb = len(self.bspace)
        for i in range(lenb):
            bspc = self.bspace[i]
            bsp1 = bspc[0]
            bsp2 = bspc[2]
            pnt1 = pnta+bsp1*vecr
            pnt2 = pnta+bsp2*vecr
            crd1 = crda+bsp1*crdr
            crd2 = crda+bsp2*crdr
            strp = LatticeStrip(lsid, pnt1, pnt2, crd1, crd2, bspc)
            ang1 = anga+bsp1*angr
            ang2 = anga+bsp2*angr
            strp.set_angles(ang1, ang2)
            strp.sht = self
            self.strps.append(strp)
            lsid += 1
        return lsid
    def set_control_points_and_normals(self):
        for pnl in self.pnls:
            pntc = pnl.return_panel_point()
            # pnl.set_coordinate(self.cord)
            pnl.set_control_point(pntc)
            pnl.set_camber_angle(0.0)
    def __repr__(self):
        return '<LatticeSheet>'

def inherit_spacing(sht: LatticeSheet):
    if sht.mirror:
        if sht.sect2.bspace is None:
            bspace = [(0.0, 0.5, 1.0)]
        else:
            lenb = len(sht.sect2.bspace)
            bspace = []
            for i in range(lenb):
                blst = []
                for j in range(2, -1, -1):
                    blst.append(1.0-sht.sect2.bspace[i][j])
                bspace.append(tuple(blst))
            bspace.reverse()
    else:
        if sht.sect1.bspace is None:
            bspace = [(0.0, 0.5, 1.0)]
        else:
            bspace = sht.sect1.bspace
    # if sht.sect1.bspace is not None:
    #     if sht.sect1.mirror:
    #         if sht.sect2.bspace is not None:
    #             bspace = [1.0-bd for bd in sht.sect2.bspace]
    #             bspace.reverse()
    #         else:
    #             bspace = [0.0, 1.0]
    #     else:
    #         bspace = sht.sect1.bspace
    # elif sht.sect2.bspace is not None:
    #     bspace = [1.0-bd for bd in sht.sect2.bspace]
    #     bspace.reverse()
    # else:
    #     bspace = [0.0, 1.0]
    # if sht.sect1.yspace is not None:
    #     if sht.sect1.mirror:
    #         if sht.sect2.yspace is None:
    #             yspace = [1.0-yd for yd in sht.sect2.yspace]
    #             yspace.reverse()
    #         else:
    #             yspace = [0.5]
    #     else:
    #         yspace = sht.sect1.yspace
    # elif sht.sect2.yspace is not None:
    #     yspace = [1.0-yd for yd in sht.sect2.yspace]
    #     yspace.reverse()
    # else:
    #     yspace = [0.5]
    return bspace
