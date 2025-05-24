from typing import TYPE_CHECKING, Any

from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d.axes3d import Axes3D
from numexpr import evaluate
from numpy import arccos, asarray, concatenate, pi, ptp, sqrt, zeros
from pygeom.geom1d import CubicSpline1D, LinearSpline1D
from pygeom.geom3d import Vector
from pygeom.tools.spacing import (equal_spacing, full_cosine_spacing,
                                  normalise_spacing)

from ..tools.airfoil import airfoil_interpolation, Airfoil
from .latticegrid import LatticeGrid
from .latticepanel import LatticePanel
from .latticesection import latticesection_from_dict
from .latticesheet import LatticeSheet
from ..tools.camber import FlatPlate

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .latticesection import LatticeSection
    from .latticestrip import LatticeStrip

    AirfoilLike = Airfoil | FlatPlate


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
    funcs: dict[str, 'SurfaceFunction'] = None

    def __init__(self, name: str, scts: list['LatticeSection'], mirror: bool,
                 funcs: dict[str, 'SurfaceFunction']) -> None:
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
        bmax = max([sct.bpos for sct in self.scts])
        for var, func in self.funcs.items():
            if var == 'twist':
                varm = '_ang'
                var1 = 'ang1'
                var2 = 'ang2'
            func.bmax = bmax
            if self.mirror:
                for i in range(hlfstrp):
                    strp = self.strps[numstrp-1-i]
                    mstrp = self.strps[i]
                    bposm = strp.bpos
                    valm = func(bposm)
                    strp.__dict__[varm] = valm
                    mstrp.__dict__[varm] = valm
                    bpos1 = strp.bpos1
                    bpos2 = strp.bpos2
                    val1 = func(bpos1)
                    val2 = func(bpos2)
                    strp.__dict__[var1] = val1
                    strp.__dict__[var2] = val2
                    mstrp.__dict__[var1] = val2
                    mstrp.__dict__[var2] = val1
            else:
                for strp in self.strps:
                    bposm = strp.bpos
                    valm = func(bposm)
                    strp.__dict__[varm] = valm
                    bpos1 = strp.bpos1
                    bpos2 = strp.bpos2
                    val1 = func(bpos1)
                    val2 = func(bpos2)
                    strp.__dict__[var1] = val1
                    strp.__dict__[var2] = val2
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

    def plot_twist(self, ax: Axes3D = None) -> Axes3D:
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        bpos = [strp.bpos for strp in self.strps]
        twist = [strp.twist for strp in self.strps]
        ax.plot(bpos, twist, label=self.name)
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


def linear_interpolate_airfoil(x: list[float],
                               af: list['AirfoilLike | None']) -> list['AirfoilLike']:
    newaf = []
    for i, afi in enumerate(af):
        if afi is None:
            a = None
            for j in range(i, -1, -1):
                if af[j] is not None:
                    a = j
                    break
            b = None
            for j in range(i, len(af)):
                if af[j] is not None:
                    b = j
                    break
            if a is None:
                xa = x[0]
                afa = FlatPlate()
            else:
                xa = x[a]
                afa = af[a]
            if b is None:
                xb = x[-1]
                afb = FlatPlate()
            else:
                xb = x[b]
                afb = af[b]
            if isinstance(afa, FlatPlate) and isinstance(afb, FlatPlate):
                afi = FlatPlate()
            elif isinstance(afa, Airfoil) and isinstance(afb, Airfoil):
                fac = (x[i] - xa)/(xb - xa)
                afi = airfoil_interpolation(afa, afb, fac)
            else:
                raise ValueError('Cannot interpolate airfoil.')
        newaf.append(afi)
    return newaf


def latticesurface_from_dict(surfdata: dict[str, Any],
                             display: bool = False) -> LatticeSurface:
    name = surfdata['name']
    mirror = surfdata.get('mirror', False)
    if display: print(f'Loading Surface: {name:s}')
    # Read in Defaults
    defaults: dict[str, Any] = surfdata.get('defaults', {})
    # Read in Functions
    funcdatas: dict[str, Any] = surfdata.get('functions', {})
    funcs = {}
    for variable, funcdata in funcdatas.items():
        funcs[variable] = surffunc_from_json(variable, funcdata)
    # Set defaults for functions to zero
    for var in funcs:
        if var not in defaults:
            defaults[var] = 0.0
    # Read Section Variables
    scts: list['LatticeSection'] = []
    for sectdata in surfdata['sections']:
        sct = latticesection_from_dict(sectdata, defaults=defaults)
        scts.append(sct)
    # Linear Interpolate Missing Variables
    x, y, z = [], [], []
    c, a, af = [], [], []
    cmb, xoc, zoc = [], [], []
    b = []
    for sct in scts:
        x.append(sct.pnt.x)
        y.append(sct.pnt.y)
        z.append(sct.pnt.z)
        c.append(sct.chord)
        a.append(sct.twist)
        af.append(sct.airfoil)
        if sct.airfoil is None:
            cmb.append(None)
        else:
            cmb.append(sct.camber)
        b.append(sct.bpos)
        xoc.append(sct.xoc)
        zoc.append(sct.zoc)
    # Check for None values in the first and last sections
    if y[0] is None or z[0] is None:
        raise ValueError('Need at least ypos or zpos specified in the first section.')
    if y[-1] is None or z[-1] is None:
        raise ValueError('Need at least ypos or zpos specified in the last section.')
    # Check for None values in the middle sections
    checky = True
    checkz = True
    for yi, zi, bi in zip(y, z, b):
        if yi is None and bi is None:
            checky = False
        if zi is None and bi is None:
            checkz = False
    # Interpolate None values in y and z
    if checky:
        linear_interpolate_none(y, z)
    elif checkz:
        linear_interpolate_none(z, y)
    else:
        raise ValueError('Need at least ypos or zpos or bpos specified in sections.')
    # Determine b values from known y and z values
    bcur = 0.0
    ycur = y[0]
    zcur = z[0]
    for i in range(len(b)):
        if b[i] is None:
            ydel = y[i] - ycur
            zdel = z[i] - zcur
            bdel = (ydel**2 + zdel**2)**0.5
            b[i] = bcur + bdel
            bcur = b[i]
            ycur = y[i]
            zcur = z[i]
    # Interpolate None values in x, y, z, c, a, cmb, xoc, zoc
    x = linear_interpolate_none(b, x)
    y = linear_interpolate_none(b, y)
    z = linear_interpolate_none(b, z)
    c = linear_interpolate_none(b, c)
    c = fill_none(c, 1.0)
    a = linear_interpolate_none(b, a)
    a = fill_none(a, 0.0)
    cmb = linear_interpolate_airfoil(b, cmb)
    xoc = linear_interpolate_none(b, xoc)
    xoc = fill_none(xoc, 0.25)
    zoc = linear_interpolate_none(b, zoc)
    zoc = fill_none(zoc, 0.0)
    display = False
    if display:
        print(f'{x = }')
        print(f'{y = }')
        print(f'{z = }')
        print(f'{c = }')
        print(f'{a = }')
        print(f'{cmb = }')
        print(f'{xoc = }')
        print(f'{zoc = }')
    for i, sct in enumerate(scts):
        sct.pnt.x = x[i]
        sct.pnt.y = y[i]
        sct.pnt.z = z[i]
        sct.chord = c[i]
        sct.twist = a[i]
        sct.camber = cmb[i]
        sct.airfoil = af[i]
        sct.xoc = xoc[i]
        sct.zoc = zoc[i]
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
    for i, (xi, yi) in enumerate(zip(x, y)):
        if yi is None and xi is None:
            continue
        elif yi is None:
            a = None
            for j in range(i, -1, -1):
                if y[j] is not None:
                    a = j
                    break
            b = None
            for j in range(i, len(y)):
                if y[j] is not None:
                    b = j
                    break
            if a is None or b is None:
                y[i] = None
            else:
                xa, xb = x[a], x[b]
                ya, yb = y[a], y[b]
                y[i] = (yb - ya)/(xb - xa)*(x[i] - xa)+ya
    return y

def fill_none(x: list[float], xval: float) -> list[float]:
    for i, xi in enumerate(x):
        if xi is None:
            x[i] = xval
    return x


class SurfaceFunction():
    variable: str | None = None
    functype: str | None = None
    bmax: float | None = None
    spline: LinearSpline1D | CubicSpline1D | None = None
    expression: str | None = None

    def __init__(self, variable: str, functype: str) -> None:
        self.variable = variable
        self.functype = functype

    def set_spline(self, spacing: str, values: list[float],
                   interp: str = 'linear') -> None:
        values: 'NDArray' = asarray(values)
        if spacing == 'equal':
            num = values.size
            nspc = equal_spacing(num - 1)
        else:
            raise ValueError('Spacing not implemented.')
        if interp == 'linear':
            self.spline = LinearSpline1D(nspc, values)
        elif interp == 'cubic':
            self.spline = CubicSpline1D(nspc, values)
        else:
            raise ValueError('Interpolation not implemented.')

    def set_expression(self, expression: str) -> None:
        self.expression = expression

    def __call__(self, value: float) -> float:
        b = value
        s = b/self.bmax
        if self.functype == 'spline':
            return self.spline.evaluate_points_at_t(s)
        elif self.functype == 'expression':
            th = arccos(s)
            local_dict = {'b': b, 's': s, 'th': th}
            global_dict = {'pi': pi}
            return evaluate(self.expression, local_dict=local_dict,
                            global_dict=global_dict)
        else:
            raise ValueError('Function type not implemented.')


def surffunc_from_json(variable: str, funcdata: dict[str, Any]) -> SurfaceFunction:
    functype = funcdata.get('functype')
    srfcfunc = SurfaceFunction(variable, functype)
    if functype == 'spline':
        splinedata: dict[str, Any] = funcdata.get('spline')
        spacing = splinedata.get('spacing', 'equal')
        interp = splinedata.get('interp', 'linear')
        values = splinedata.get('values')
        srfcfunc.set_spline(spacing, values, interp)
    elif functype == 'expression':
        expression = funcdata.get('expression')
        srfcfunc.set_expression(expression)
    else:
        raise ValueError('Function type not implemented.')
    return srfcfunc
