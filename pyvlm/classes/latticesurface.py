from math import sqrt
from pygeom.matrix3d import zero_matrix_vector
from matplotlib.pyplot import figure
from .latticesheet import LatticeSheet
from .latticepanel import LatticePanel

class LatticeSurface(object):
    name = None
    scts = None
    shts = None
    cspc = None
    xspace = None
    strps = None
    pnts = None
    pnls = None
    area = None
    sgrp = None
    def __init__(self, name: str, scts: list, mirror: bool, funcs: list):
        self.name = name
        self.scts = scts
        self.mirror = mirror
        self.funcs = funcs
        self.update()
    def update(self):
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
    def set_chord_distribution(self, cspc: list):
        from pyvlm.tools import normalise_spacing
        self.cspc = normalise_spacing(cspc)
    def set_chord_equal_distribution(self, cnum: int):
        from pyvlm.tools import equal_spacing
        csp = equal_spacing(4*cnum)
        self.cspc = [tuple(csp[i*4:i*4+5]) for i in range(cnum)]
    def set_chord_cosine_distribution(self, cnum: int):
        from pyvlm.tools import full_cosine_spacing
        if cnum > 1:
            csp = full_cosine_spacing(4*cnum+2)
            csp = [0.0]+csp[2:-2]+[1.0]
            self.cspc = [tuple(csp[i*4:i*4+5]) for i in range(cnum)]
        else:
            self.set_chord_equal_distribution(cnum)
    def mesh(self, lsid: int, lpid: int):
        from pygeom.geom3d import Vector
        from numpy.matlib import empty
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
        self.pnts = empty((lenb, lenc+1), dtype=Vector)
        for i in range(lenb):
            minx = pnts[i].x
            y = pnts[i].y
            z = pnts[i].z
            c = crds[i]
            cd = self.cspc[0][0]
            x = minx+cd*c
            self.pnts[i, 0] = Vector(x, y, z)
            for j in range(1, lenc+1):
                cd = self.cspc[j-1][-1]
                x = minx+cd*c
                self.pnts[i, j] = Vector(x, y, z)
        self.pnls = empty((lenb-1, lenc), dtype=LatticePanel)
        for i, strp in enumerate(self.strps):
            for j in range(lenc):
                pnts = [
                    self.pnts[i, j],
                    self.pnts[i+1, j],
                    self.pnts[i, j+1],
                    self.pnts[i+1, j+1]
                ]
                cspc = self.cspc[j]
                pnl = LatticePanel(lpid, pnts, cspc, strp)
                self.pnls[i, j] = pnl
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
    def plot_surface(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca(projection='3d')
            ax.grid(True)
        x, y, z = self.point_xyz()
        ax.plot_surface(x, y, z, label=self.name)
        return ax
    @property
    def strpb(self):
        return [strp.bpos for strp in self.strps]
    @property
    def strpy(self):
        return [strp.pnti.y for strp in self.strps]
    @property
    def strpz(self):
        return [strp.pnti.z for strp in self.strps]
    @property
    def strpi(self):
        return [strp.lsid for strp in self.strps]
    @property
    def lstrpi(self):
        return self.sgrp[0]
    @property
    def mstrpi(self):
        return self.sgrp[1]
    @property
    def pnli(self):
        lpids = []
        for i in range(self.pnls.shape[0]):
            for j in range(self.pnls.shape[1]):
                lpids.append(self.pnls[i, j].lpid)
        return lpids
    def vortex_line_points(self, indp: int, nump: int):
        nums = len(self.strps)
        num = nums*nump+1
        rpt = zero_matrix_vector((num, 1))
        j = 0
        for strp in self.strps:
            pnl = strp.pnls[indp]
            for i in range(nump):
                pnt = pnl.pnta+i/nump*pnl.leni
                rpt[j, 0] = pnt
                j += 1
        rpt[j, 0] = pnl.pntb
        return rpt
    def __repr__(self):
        return '<LatticeSurface {:s}>'.format(self.name)

def latticesurface_from_json(surfdata: dict, display: bool=False):
    from .latticesection import latticesection_from_json
    name = surfdata['name']
    if 'mirror' in surfdata:
        mirror = surfdata['mirror']
    else:
        mirror = False
    if display: print(f'Loading Surface: {name:s}')
    # Read Section Variables
    scts = []
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
            return ValueError
        else:
            y = linear_interpolate_none(z, y)
    else:
        z = linear_interpolate_none(y, z)
    lenscts = len(scts)
    b = [0.0]
    for i in range(lenscts-1):
        bi = b[i]+sqrt((y[i+1]-y[i])**2+(z[i+1]-z[i])**2)
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
    xpos, ypos, zpos = 0.0, 0.0, 0.0
    if 'xpos' in surfdata:
        xpos = surfdata['xpos']
    if 'ypos' in surfdata:
        ypos = surfdata['ypos']
    if 'zpos' in surfdata:
        zpos = surfdata['zpos']
    twist = 0.0
    if 'twist' in surfdata:
        twist = surfdata['twist']
    if 'ruled' in surfdata:
        ruled = surfdata['ruled']
    else:
        ruled = False
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

def linear_interpolate_none(x: list, y: list):
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

class SurfaceFunction(object):
    var = None
    dist = None
    interp = None
    values = None
    spline = None
    def __init__(self, var: str, spacing: str, interp: str, values: list):
        self.var = var
        self.spacing = spacing
        self.interp = interp
        self.values = values
    def set_spline(self, bmax: float):
        if self.spacing == 'equal':
            num = len(self.values)
            from pyvlm.tools import equal_spacing
            nspc = equal_spacing(num-1)
            spc = [bmax*nspci for nspci in nspc]
        if self.interp == 'linear':
            from pygeom.geom1d import LinearSpline
            self.spline = LinearSpline(spc, self.values)
        elif self.interp == 'cubic':
            from pygeom.geom1d import CubicSpline
            self.spline = CubicSpline(spc, self.values)
    def interpolate(self, b: float):
        return self.spline.single_interpolate_spline(b)

def surffunc_from_json(funcdata: dict):
    var = funcdata["variable"]
    if "spacing" in funcdata:
        spacing = funcdata["spacing"]
    else:
        spacing = "equal"
    if "interp" in funcdata:
        interp = funcdata["interp"]
    else:
        interp = "linear"
    values = funcdata["values"]
    return SurfaceFunction(var, spacing, interp, values)
