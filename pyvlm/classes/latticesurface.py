from math import pi, cos, sqrt
from .latticestrip import LatticeStrip
from .latticesheet import LatticeSheet
from .latticepanel import LatticePanel
from pygeom.matrix3d import zero_matrix_vector
from matplotlib.pyplot import figure

class LatticeSurface(object):
    name = None
    sects = None
    cspace = None
    xspace = None
    strps = None
    pnts = None
    pnls = None
    msect = None
    area = None
    def __init__(self, name: str, sects: list, mirror: bool, funcs:list):
        self.name = name
        self.sects = sects
        self.mirror = mirror
        self.funcs = funcs
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
        csp = equal_spacing(4*numc)
        self.cspace = [tuple(csp[i*4:i*4+5]) for i in range(numc)]
    def set_chord_cosine_distribution(self, numc: int):
        from pyvlm.tools import full_cosine_spacing
        if numc > 1:
            csp = full_cosine_spacing(4*numc+2)
            csp = [0.0]+csp[2:-2]+[1.0]
            self.cspace = [tuple(csp[i*4:i*4+5]) for i in range(numc)]
        else:
            self.set_chord_equal_distribution(numc)
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
        self.pnts = empty((lenb, lenc+1), dtype=Point)
        for i in range(lenb):
            minx = pnts[i].x
            y = pnts[i].y
            z = pnts[i].z
            c = crds[i]
            cd = self.cspace[0][0]
            x = minx+cd*c
            self.pnts[i, 0] = Point(x, y, z)
            for j in range(1, lenc+1):
                cd = self.cspace[j-1][-1]
                x = minx+cd*c
                self.pnts[i, j] = Point(x, y, z)
        self.pnls = empty((lenb-1, lenc), dtype=LatticePanel)
        for i, strp in enumerate(self.strps):
            for j in range(lenc):
                pnts = [
                    self.pnts[i, j],
                    self.pnts[i+1, j],
                    self.pnts[i, j+1],
                    self.pnts[i+1, j+1]
                ]
                cspc = self.cspace[j]
                pnl = LatticePanel(lpid, pnts, cspc, strp)
                self.pnls[i, j] = pnl
                lpid += 1
        if self.mirror:
            numstrp = len(self.strps)
            hlfstrp = int(numstrp/2)
            for i in range(hlfstrp):
                strp = self.strps[numstrp-1-i]
                mstrp = self.strps[i]
                strp.set_mirror(mstrp)
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
        for i, sect in enumerate(self.sects):
            sect.bpos = bpos[i]
        for sht in self.shts:
            sht.set_strip_bpos()
        bmax = max(bpos)
        for func in self.funcs:
            func.set_spline(bmax)
            var = func.var
            if var == 'angle':
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
        return [strp.lsid for strp in self.strps if strp.msid is None]
    @property
    def mstrpi(self):
        return [strp.lsid for strp in self.strps if strp.msid is not None]
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
    from .latticesection import latticesecttion_from_json
    name = surfdata['name']
    if 'mirror' in surfdata:
        mirror = surfdata['mirror']
    else:
        mirror = False
    if display: print(f'Loading Surface: {name:s}')
    # Read Section Variables
    sects = []
    for sectdata in surfdata['sections']:
        sect = latticesecttion_from_json(sectdata)
        sects.append(sect)
    # Linear Interpolate Missing Variables
    x, y, z, c, a = [], [], [], [], []
    for sect in sects:
        x.append(sect.pnt.x)
        y.append(sect.pnt.y)
        z.append(sect.pnt.z)
        c.append(sect.chord)
        a.append(sect.angle)
    if None in y:
        if None is z:
            return ValueError
        else:
            y = linear_interpolate_none(z, y)
    else:
        z = linear_interpolate_none(y, z)
    lensects = len(sects)
    b = [0.0]
    for i in range(lensects-1):
        bi = b[i]+sqrt((y[i+1]-y[i])**2+(z[i+1]-z[i])**2)
        b.append(bi)
    x = linear_interpolate_none(b, x)
    c = linear_interpolate_none(b, c)
    a = linear_interpolate_none(b, a)
    for i, sect in enumerate(sects):
        sect.pnt.x = x[i]
        sect.pnt.y = y[i]
        sect.pnt.z = z[i]
        sect.chord = c[i]
        sect.angle = a[i]
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
    angle = 0.0
    if 'angle' in surfdata:
        angle = surfdata['angle']
    for sect in sects:
        sect.offset_position(xpos, ypos, zpos)
        sect.offset_angle(angle)
    surf = LatticeSurface(name, sects, mirror, funcs)
    if 'numc' in surfdata and 'cspace' in surfdata:
        numc = surfdata['numc']
        cspace = surfdata['cspace']
        if cspace == 'equal':
            surf.set_chord_equal_distribution(numc)
        elif cspace == 'cosine':
            surf.set_chord_cosine_distribution(numc)
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
