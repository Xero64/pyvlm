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
    ctrls = None
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
        self.ctrls = inherit_controls(self)
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
    def inherit_panels(self):
        self.pnls = []
        for strp in self.strps:
            for pnl in strp.pnls:
                self.pnls.append(pnl)
    def set_control_panels(self):
        # print(f'# Sheet Panels = {len(self.pnls):d}')
        for control in self.ctrls:
            # print(f'Setting {control:s} panels!')
            ctrl = self.ctrls[control]
            for pnl in self.pnls:
                if pnl.cspc[3] >= ctrl.xhinge:
                    ctrl.add_panel(pnl)
    # def set_control_points_and_normals(self):
    #     for pnl in self.pnls:
    #         pntc = pnl.return_panel_point()
    #         pnl.set_control_point(pntc)
    #         pnl.set_camber_angle(0.0)
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
    return bspace

def inherit_controls(sht: LatticeSheet):
    ctrls = {}
    if sht.mirror:
        for control in sht.sect2.ctrls:
            ctrl = sht.sect2.ctrls[control]
            newctrl = ctrl.duplicate(mirror=True)
            ctrls[control] = newctrl
    else:
        for control in sht.sect1.ctrls:
            ctrl = sht.sect1.ctrls[control]
            newctrl = ctrl.duplicate(mirror=False)
            ctrls[control] = newctrl
    for control in ctrls:
        ctrl = ctrls[control]
        if ctrl.uhvec.return_magnitude() == 0.0:
            pnt1 = sht.sect1.pnt
            crd1 = sht.sect1.chord
            pnta = pnt1+crd1*ihat*ctrl.xhinge
            pnt2 = sht.sect2.pnt
            crd2 = sht.sect2.chord
            pntb = pnt2+crd2*ihat*ctrl.xhinge
            hvec = pntb-pnta
            ctrl.set_hinge_vector(hvec)
    return ctrls
