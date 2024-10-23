from typing import TYPE_CHECKING

from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d.axes3d import Axes3D
from numpy import asarray, concatenate, ptp, sqrt, zeros
from pygeom.geom1d import CubicSpline1D, LinearSpline1D
from pygeom.geom3d import Vector
from pygeom.tools.spacing import (equal_spacing, full_cosine_spacing,
                                  normalise_spacing)

from .latticegrid import LatticeGrid
from .latticepanel import LatticePanel
from .latticesection import latticesection_from_json
from .latticesheet import LatticeSheet

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .latticesection import LatticeSection
    from .latticestrip import LatticeStrip


class LatticeSurface():
    name: str = None
    scts: list['LatticeSection'] = None
    shts: list[LatticeSheet] = None
    cspc: list[float] = None
    strps: list['LatticeStrip'] = None
    pnts: list[list[LatticeGrid]] = None
    pnls: list[list[LatticePanel]] = None
    area: float = None
    sgrp: list[list[int]] = None
    funcs: list['SurfaceFunction'] = None

    def __init__(self, name: str, scts: list, mirror: bool, funcs: list) -> None:
        self.name = name
        self.scts = scts
        self.mirror = mirror
        self.funcs = funcs
        self.update()

    def update(self) -> None:
        if self.mirror and self.scts[0].pnt.y == 0.0:
            numsct = len(self.scts)
            newscts = []
            for i in range(numsct-1):
                sct = self.scts[numsct-1-i]
                msct = sct.return_mirror()
                newscts.append(msct)
            for sct in self.scts:
                newscts.append(sct)
            self.scts = newscts
        elif self.mirror and self.scts[0].pnt.y != 0.0:
            print(f'Warning: Cannot mirror {self.name}.')
            self.mirror = False

    def set_chord_distribution(self, cspc: list[float]) -> None:
        self.cspc = normalise_spacing(cspc)

    def set_chord_equal_distribution(self, cnum: int) -> None:
        csp = equal_spacing(4*cnum)
        self.cspc = [csp[i*4:i*4 + 5] for i in range(cnum)]

    def set_chord_cosine_distribution(self, cnum: int) -> None:
        if cnum > 1:
            csp = full_cosine_spacing(4*cnum + 2)
            csp = concatenate(([0.0], csp[2:-2], [1.0]))
            self.cspc = [csp[i*4:i*4 + 5] for i in range(cnum)]
        else:
            self.set_chord_equal_distribution(cnum)

    def mesh(self, lsid: int, lpid: int) -> None:
        nums = len(self.scts)
        self.shts = []
        for i in range(nums-1):
            a, b = i, i+1
            scta = self.scts[a]
            sctb = self.scts[b]
            self.shts.append(LatticeSheet(scta, sctb))
        self.strps = []
        for sht in self.shts:
            lsid = sht.mesh_strips(lsid)
            self.strps += sht.strps
        pnts = [strp.pnt1 for strp in self.strps]
        pnts.append(self.strps[-1].pnt2)
        crds = [strp.crd1 for strp in self.strps]
        crds.append(self.strps[-1].crd2)
        lenb = len(pnts)
        lenc = len(self.cspc)
        self.pnts = []
        for i in range(lenb):
            c = crds[i]
            cd = self.cspc[0][0]
            minx = pnts[i].x
            x = minx + cd*c
            y = pnts[i].y
            z = pnts[i].z
            self.pnts.append([])
            self.pnts[i].append(LatticeGrid(x, y, z))
            for j in range(1, lenc+1):
                cd = self.cspc[j-1][-1]
                x = minx + cd*c
                self.pnts[i].append(LatticeGrid(x, y, z))
        self.pnls = []
        for i, strp in enumerate(self.strps):
            self.pnls.append([])
            for j in range(lenc):
                pnts = [
                    self.pnts[i][j],
                    self.pnts[i+1][j],
                    self.pnts[i][j+1],
                    self.pnts[i+1][j+1]
                ]
                cspc = self.cspc[j]
                pnl = LatticePanel(lpid, pnts, cspc, strp)
                self.pnls[i].append(pnl)
                lpid += 1
        if self.mirror:
            self.sgrp = [[], []]
            numstrp = len(self.strps)
            hlfstrp = int(numstrp/2)
            for i in range(hlfstrp):
                lstrp = self.strps[numstrp-1-i]
                mstrp = self.strps[i]
                self.sgrp[0].append(lstrp.lsid)
                self.sgrp[1].append(mstrp.lsid)
        else:
            self.sgrp = [[]]
            numstrp = len(self.strps)
            for i in range(numstrp):
                lstrp = self.strps[numstrp-1-i]
                self.sgrp[0].append(lstrp.lsid)
        bpos = [0.0]
        for sht in self.shts:
            sht.inherit_panels()
            sht.set_control_panels()
            bpos.append(bpos[-1]+sht.width)
        if self.mirror:
            numsht = len(self.shts)
            wmir = bpos[int(numsht/2)]
            for i in range(len(bpos)):
                bpos[i] = bpos[i]-wmir
        for i, sct in enumerate(self.scts):
            sct.bpos = bpos[i]
        for sht in self.shts:
            sht.set_strip_bpos()
        bmax = max(bpos)
        for func in self.funcs:
            func.set_spline(bmax)
            var = func.var
            if var == 'twist':
                var = '_ang'
            if self.mirror:
                for i in range(hlfstrp):
                    strp = self.strps[numstrp-1-i]
                    mstrp = self.strps[i]
                    bpos = strp.bpos
                    val = func.interpolate(bpos)
                    strp.__dict__[var] = val
                    mstrp.__dict__[var] = val
            else:
                for strp in self.strps:
                    bpos = strp.bpos
                    val = func.interpolate(bpos)
                    strp.__dict__[var] = val
        self.area = 0.0
        for sht in self.shts:
            if not sht.noload:
                self.area += sht.area
        return lsid, lpid

    def point_xyz(self) -> tuple['NDArray', 'NDArray', 'NDArray']:
        shape = (len(self.pnts), len(self.pnts[0]))
        x = zeros(shape)
        y = zeros(shape)
        z = zeros(shape)
        for i in range(len(self.pnts)):
            for j in range(len(self.pnts[i])):
                pnt = self.pnts[i][j]
                x[i, j] = pnt.x
                y[i, j] = pnt.y
                z[i, j] = pnt.z
        return x, y, z

    def return_panels(self) -> list[LatticePanel]:
        pnls = []
        for i in range(len(self.pnls)):
            for j in range(len(self.pnls[i])):
                pnls.append(self.pnls[i][j])
        return pnls

    def plot_surface(self, ax: Axes3D = None) -> Axes3D:
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = Axes3D(fig)
            fig.add_axes(ax)
            ax.grid(True)
        x, y, z = self.point_xyz()
        ax.set_box_aspect((ptp(x), ptp(y), ptp(z)))
        ax.plot_surface(x, y, z, label=self.name)
        return ax

    @property
    def strpb(self) -> list['LatticeStrip']:
        return [strp.bpos for strp in self.strps]

    @property
    def strpy(self) -> list['LatticeStrip']:
        return [strp.pnti.y for strp in self.strps]

    @property
    def strpz(self) -> list['LatticeStrip']:
        return [strp.pnti.z for strp in self.strps]

    @property
    def strpi(self) -> list['LatticeStrip']:
        return [strp.lsid for strp in self.strps]

    @property
    def lstrpi(self) -> list[int]:
        return self.sgrp[0]

    @property
    def mstrpi(self) -> list[int]:
        return self.sgrp[1]

    @property
    def pnli(self) -> list[int]:
        lpids = []
        for i in range(len(self.pnls)):
            for j in range(len(self.pnls[i])):
                lpids.append(self.pnls[i][j].lpid)
        return lpids

    def vortex_line_points(self, indp: int, nump: int) -> Vector:
        nums = len(self.strps)
        num = nums*nump+1
        rpt = Vector.zeros((num, 1))
        j = 0
        for strp in self.strps:
            pnl = strp.pnls[indp]
            for i in range(nump):
                pnt = pnl.pnta+i/nump*pnl.leni
                rpt[j, 0] = pnt
                j += 1
        rpt[j, 0] = pnl.pntb
        return rpt

    def __repr__(self) -> str:
        return '<LatticeSurface {:s}>'.format(self.name)


def latticesurface_from_json(surfdata: dict, display: bool=False) -> LatticeSurface:
    name = surfdata['name']
    mirror = surfdata.get('mirror', False)
    if display: print(f'Loading Surface: {name:s}')
    # Read Section Variables
    scts: list['LatticeSection'] = []
    for sectdata in surfdata['sections']:
        sct = latticesection_from_json(sectdata)
        scts.append(sct)
    # Linear Interpolate Missing Variables
    x, y, z, c, a = [], [], [], [], []
    for sct in scts:
        x.append(sct.pnt.x)
        y.append(sct.pnt.y)
        z.append(sct.pnt.z)
        c.append(sct.chord)
        a.append(sct.twist)
    if None in y:
        if None is z:
            raise ValueError('Either y or z must be defined.')
        else:
            y = linear_interpolate_none(z, y)
    else:
        z = linear_interpolate_none(y, z)
    lenscts = len(scts)
    b = [0.0]
    for i in range(lenscts-1):
        bi = b[i] + sqrt((y[i+1]-y[i])**2+(z[i+1]-z[i])**2)
        b.append(bi)
    x = linear_interpolate_none(b, x)
    c = linear_interpolate_none(b, c)
    a = linear_interpolate_none(b, a)
    for i, sct in enumerate(scts):
        sct.pnt.x = x[i]
        sct.pnt.y = y[i]
        sct.pnt.z = z[i]
        sct.chord = c[i]
        sct.twist = a[i]
    # Read in Function Data
    funcs = []
    if 'functions' in surfdata:
        for funcdata in surfdata['functions']:
            func = surffunc_from_json(funcdata)
            funcs.append(func)
    # Entire Surface Position
    xpos = surfdata.get('xpos', 0.0)
    ypos = surfdata.get('ypos', 0.0)
    zpos = surfdata.get('zpos', 0.0)
    twist = surfdata.get('twist', 0.0)
    ruled = surfdata.get('ruled', False)
    for sct in scts:
        sct.offset_position(xpos, ypos, zpos)
        sct.offset_twist(twist)
        sct.ruled = ruled
    surf = LatticeSurface(name, scts, mirror, funcs)
    if 'cnum' in surfdata:
        cnum = surfdata['cnum']
        cspc = 'cosine'
        if 'cspc' in surfdata:
            cspc = surfdata['cspc'].lower()
        if cspc == 'equal':
            surf.set_chord_equal_distribution(cnum)
        elif cspc in ('cosine', 'full-cosine'):
            surf.set_chord_cosine_distribution(cnum)
    return surf

def linear_interpolate_none(x: list[float], y: list[float]) -> list[float]:
    for i, yi in enumerate(y):
        if yi is None:
            for j in range(i, -1, -1):
                if y[j] is not None:
                    a = j
                    break
            for j in range(i, len(y)):
                if y[j] is not None:
                    b = j
                    break
            xa, xb = x[a], x[b]
            ya, yb = y[a], y[b]
            y[i] = (yb-ya)/(xb-xa)*(x[i]-xa)+ya
    return y


class SurfaceFunction():
    var: str = None
    interp: str = None
    values: 'NDArray' = None
    spline: LinearSpline1D | CubicSpline1D = None

    def __init__(self, var: str, spacing: str, interp: str,
                 values: list[float]) -> None:
        self.var = var
        self.spacing = spacing
        self.interp = interp
        self.values = asarray(values)

    def set_spline(self, bmax: float) -> None:
        if self.spacing == 'equal':
            num = len(self.values)
            nspc = equal_spacing(num-1)
        if self.interp == 'linear':
            self.spline = LinearSpline1D(nspc, self.values)
        elif self.interp == 'cubic':
            self.spline = CubicSpline1D(nspc, self.values)

    def interpolate(self, value: float) -> float:
        return self.spline.evaluate_points_at_t(value)

def surffunc_from_json(funcdata: dict) -> SurfaceFunction:
    variable = funcdata.get('variable')
    spacing = funcdata.get('spacing', 'equal')
    interp = funcdata.get('interp', 'linear')
    values = funcdata.get('values')
    return SurfaceFunction(variable, spacing, interp, values)
