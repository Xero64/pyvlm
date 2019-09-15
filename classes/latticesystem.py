from math import pi, atan2
from numpy.matlib import zeros, empty
from numpy.linalg import solve, inv
from pygeom.geom3d import Point, Vector, ihat, jhat, zero_vector

class LatticeSystem(object):
    name = None
    srfcs = None
    strps = None
    mstrp = None
    pnls = None
    mpnl = None
    bref = None
    cref = None
    sref = None
    rref = None
    _gam = None
    _avg = None
    _afg = None
    _amg = None
    _avc = None
    _aic = None
    _afs = None
    _adc = None
    _bvg = None
    _bdg = None
    _blg = None
    _byg = None
    _bmg = None
    _ar = None
    _strpy = None
    def __init__(self, name: str, srfcs: list):
        self.name = name
        self.srfcs = srfcs
        self.mesh()
    def mesh(self):
        lsid = 0
        lpid = 0
        for srfc in self.srfcs:
            lsid, lpid = srfc.mesh(lsid, lpid)
        self.inherit()
    def inherit(self):
        strpdct = {}
        pnldct = {}
        istrp = 0
        ipnl = 0
        self.strps = []
        self.mstrp = []
        self.pnls = []
        self.mpnl = []
        for srfc in self.srfcs:
            for strp in srfc.strps:
                self.strps.append(strp)
                if strp.msid is None:
                    strpind = istrp
                    strpdct[strp.lsid] = istrp
                else:
                    strpind = strpdct[strp.msid]
                self.mstrp.append(strpind)
                istrp += 1
            pnls = srfc.return_panels()
            for pnl in pnls:
                self.pnls.append(pnl)
                if pnl.mpid is None:
                    pnlind = ipnl
                    pnldct[pnl.lpid] = ipnl
                else:
                    pnlind = pnldct[pnl.mpid]
                self.mpnl.append(pnlind)
                ipnl += 1
    def reset(self):
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None
    def set_reference_geometry(self, bref: float, cref: float, sref: float):
        self.bref = bref
        self.cref = cref
        self.sref = sref
        self._ar = None
    def set_reference_point(self, xref: float, yref: float, zref: float):
        self.rref = Point(xref, yref, zref)
    @property
    def gam(self):
        if self._gam is None:
            self._gam = -solve(self.aic, self.afs)
        return self._gam
    @property
    def avg(self):
        if self._avg is None:
            num = len(self.pnls)
            self._avg = empty((num, num), dtype=Vector)
            for i, pnli in enumerate(self.pnls):
                for j, pnlj in enumerate(self.pnls):
                    self._avg[i, j] = pnlj.velocity(pnli.pnti)
        return self._avg
    @property
    def afg(self):
        if self._afg is None:
            num = len(self.pnls)
            self._afg = empty((num, num), dtype=Vector)
            for i, pnl in enumerate(self.pnls):
                for j in range(num):
                    self._afg[i, j] = self.avg[i, j]**pnl.leni
        return self._afg
    @property
    def amg(self):
        if self._amg is None:
            num = len(self.pnls)
            self._amg = empty((num, num), dtype=Vector)
            for i, pnl in enumerate(self.pnls):
                for j in range(num):
                    self._amg[i, j] = (pnl.pnti-self.rref)**self.afg[i, j]
        return self._amg
    @property
    def avc(self):
        if self._avc is None:
            num = len(self.pnls)
            self._avc = empty((num, num), dtype=Vector)
            for i, pnli in enumerate(self.pnls):
                for j, pnlj in enumerate(self.pnls):
                    self._avc[i, j] = pnlj.velocity(pnli.pntc)
        return self._avc
    @property
    def aic(self):
        if self._aic is None:
            num = len(self.pnls)
            self._aic = zeros((num, num))
            for i, pnl in enumerate(self.pnls):
                for j in range(num):
                    self._aic[i, j] = self.avc[i, j]*pnl.nrml
        return self._aic
    @property
    def adc(self):
        if self._adc is None:
            num = len(self.pnls)
            self._adc = zeros((num, num))
            for i, pnl in enumerate(self.pnls):
                for j in range(num):
                    self._adc[i, j] = self.avc[i, j]*pnl.tang
        return self._adc
    @property
    def afs(self):
        if self._afs is None:
            num = len(self.pnls)
            self._afs = zeros((num, 3))
            for i, pnl in enumerate(self.pnls):
                self._afs[i, 0] = pnl.nrml.x
                self._afs[i, 1] = pnl.nrml.y
                self._afs[i, 2] = pnl.nrml.z
        return self._afs
    @property
    def bvg(self):
        if self._bvg is None:
            num = len(self.strps)
            self._bvg = zeros((num, num))
            for i, strpi in enumerate(self.strps):
                for j, strpj in enumerate(self.strps):
                    self._bvg[i, j] = strpj.trefftz_velocity(strpi.pnti)*strpi.nrmt
        return self._bvg
    @property
    def bdg(self):
        if self._bdg is None:
            num = len(self.strps)
            self._bdg = zeros((num, num))
            for i, strp in enumerate(self.strps):
                for j in range(num):
                    self._bdg[i, j] = -strp.dst*self.bvg[i, j]/2
        return self._bdg
    @property
    def blg(self):
        if self._blg is None:
            num = len(self.strps)
            self._blg = zeros((num, 1))
            for i, strp in enumerate(self.strps):
                self._blg[i, 0] = strp.lent.y
        return self._blg
    @property
    def byg(self):
        if self._byg is None:
            num = len(self.strps)
            self._byg = zeros((num, 1))
            for i, strp in enumerate(self.strps):
                self._byg[i, 0] = -strp.lent.z
        return self._byg
    @property
    def bmg(self):
        if self._bmg is None:
            num = len(self.strps)
            self._bmg = zeros((num, 1))
            for i, strp in enumerate(self.strps):
                self._bmg[i, 0] = strp.pnti.y*self.blg[i, 0]-strp.pnti.z*self.byg[i, 0]
        return self._bmg
    @property
    def strpy(self):
        if self._strpy is None:
            self._strpy = [strp.pnti.y for strp in self.strps]
        return self._strpy
    # def optimum_lift_distribution(self, Lspec: float, lspec: float=None):
        # from IPython.display import display
        # from pyhtml import HTMLMatrix
        # nump = len(self.strps)
        # numc = 1
        # if lspec is not None:
        #     numc += 2
        # amat = zeros((nump+numc, nump+numc))
        # bmat = zeros((nump+numc, 1))
        # amat[0:nump, 0:nump] = self.bdg + self.bdg.transpose()
        # amat[0:nump, nump] = self.blg
        # amat[nump, 0:nump] = self.blg.transpose()
        # bmat[nump, 0] = Lspec
        # bmg1 = self.bmg.copy()
        # bmg2 = self.bmg.copy()
        # hlfp = int(nump/2)
        # for i in range(hlfp):
        #     bmg1[i, 0] = 0.0
        # for i in range(hlfp, nump):
        #     bmg2[i, 0] = 0.0
        # if lspec is not None:
        #     amat[0:nump, nump+1] = bmg1
        #     amat[nump+1, 0:nump] = bmg1.transpose()
        #     bmat[nump+1, 0] = lspec
        #     amat[0:nump, nump+2] = bmg2
        #     amat[nump+2, 0:nump] = bmg2.transpose()
        #     bmat[nump+2, 0] = -1.0*lspec
        # xmat = solve(amat, bmat)
        # phi = [xmat[i, 0] for i in range(nump)]
        # lam = [xmat[nump+j, 0] for j in range(numc)]
        # Di = (xmat[0:nump, 0].transpose()*self.bdg*xmat[0:nump, 0])[0, 0]
        # L = (xmat[0:nump, 0].transpose()*self.blg)[0, 0]
        # l = (xmat[0:nump, 0].transpose()*bmg1)[0, 0]
        # amat = HTMLMatrix(amat)
        # display(amat)
        # return phi, lam, Di, L, l
    @property
    def ar(self):
        if self._ar is None:
            self._ar = self.bref**2/self.sref
        return self._ar
    def set_strip_alpha(self, alpha: list):
        for i, strp in enumerate(self.strps):
            strp._ang = alpha[i]
        self._aic = None
        self._afs = None
        self._gam = None
    def plot_surface_ipv(self):
        from ipyvolume import figure, plot_wireframe, show
        fig = figure()
        for srfc in self.srfcs:
            x, y, z = srfc.point_xyz()
            _ = plot_wireframe(x, y, z)
        show()
        return fig
    def plot_surface(self, view=None):
        from matplotlib.pyplot import figure
        from mpl_toolkits.mplot3d import axes3d
        # pltbox = PlotBox3D()
        fig = figure(figsize=(12,8))
        ax = fig.gca(projection='3d')
        ax.set_proj_type('ortho')
        ax.grid(False)
        # ax.set_aspect('equal')
        for srfc in self.srfcs:
            x, y, z = srfc.point_xyz()
            ax.plot_wireframe(x, y, z)
            # pltbox.update_box(x, y, z)
        if view == 'top':
            ax.view_init(90.0, 0.0)
        elif view == 'left':
            ax.view_init(0.0, 0.0)
        # pltbox.plot_box(ax)
        set_axes_equal(ax)
        fig.tight_layout()
        return ax
    def print_strip_geometry(self):
        from py2md.classes import MDTable
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('Xle', '.5f')
        table.add_column('Yle', '.5f')
        table.add_column('Zle', '.5f')
        table.add_column('Chord', '.4f')
        table.add_column('Width', '.5f')
        table.add_column('Area', '.6f')
        table.add_column('Dihed', '.4f')
        table.add_column('Incid', '.4f')
        for strp in self.strps:
            j = strp.lsid
            xle = strp.pnti.x
            yle = strp.pnti.y
            zle = strp.pnti.z
            chord = strp.chord
            width = strp.dst
            area = strp.area
            dihed = strp.dihedral
            angle = strp.angle
            table.add_row([j, xle, yle, zle, chord, width, area, dihed, angle])
        print(table)
    def print_panel_geometry(self):
        from py2md.classes import MDTable
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('X', '.5f')
        table.add_column('Y', '.5f')
        table.add_column('Z', '.5f')
        table.add_column('DX', '.5f')
        for pnl in self.pnls:
            j = pnl.lpid
            x = pnl.pnti.x
            y = pnl.pnti.y
            z = pnl.pnti.z
            dx = pnl.crd
            table.add_row([j, x, y, z, dx])
        print(table)
    def __repr__(self):
        return '<LatticeSystem: {:s}>'.format(self.name)
    def __str__(self):
        from py2md.classes import MDTable
        outstr = '# Lattice System '+self.name+'\n'
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('Sref', 'g', data=[self.sref])
        table.add_column('cref', 'g', data=[self.cref])
        table.add_column('bref', 'g', data=[self.bref])
        table.add_column('xref', '.3f', data=[self.rref.x])
        table.add_column('yref', '.3f', data=[self.rref.y])
        table.add_column('zref', '.3f', data=[self.rref.z])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('# Strips', 'd', data=[len(self.strps)])
        table.add_column('# Panels', 'd', data=[len(self.pnls)])
        outstr += table._repr_markdown_()
        return outstr
    def _repr_markdown_(self):
        return self.__str__()

def latticesystem_from_json(jsonfilepath: str):
    from .latticesurface import latticesurface_from_json
    from json import load
    with open(jsonfilepath, 'rt') as jsonfile:
        data = load(jsonfile)
    name = data['name']
    sfcs = []
    for surfdata in data['surfaces']:
        sfc = latticesurface_from_json(surfdata)
        sfcs.append(sfc)
    sys = LatticeSystem(name, sfcs)
    bref = data['bref']
    cref = data['cref']
    sref = data['sref']
    sys.set_reference_geometry(bref, cref, sref)
    xref = data['xref']
    yref = data['yref']
    zref = data['zref']
    sys.set_reference_point(xref, yref, zref)
    return sys

class PlotBox3D(object):
    xmin = None
    ymin = None
    zmin = None
    xmax = None
    ymax = None
    zmax = None
    def __init__(self):
        pass
    def update_box(self, x, y, z):
        if self.xmin is None:
            self.xmin = x.min()
        else:
            self.xmin = min(x.min(), self.xmin)
        if self.xmax is None:
            self.xmax = x.max()
        else:
            self.xmax = max(x.max(), self.xmax)
        if self.ymin is None:
            self.ymin = y.min()
        else:
            self.ymin = min(y.min(), self.ymin)
        if self.ymax is None:
            self.ymax = y.max()
        else:
            self.ymax = max(y.max(), self.ymax)
        if self.zmin is None:
            self.zmin = z.min()
        else:
            self.zmin = min(z.min(), self.zmin)
        if self.zmax is None:
            self.zmax = z.max()
        else:
            self.zmax = max(z.max(), self.zmax)
    def plot_box(self, ax):
        xrng = self.xmax-self.xmin
        yrng = self.ymax-self.ymin
        zrng = self.zmax-self.zmin
        xctr = (self.xmax+self.xmin)/2
        yctr = (self.ymax+self.ymin)/2
        zctr = (self.zmax+self.zmin)/2
        maxrng = max([xrng, yrng, zrng])
        xmin = xctr-maxrng/2
        ymin = yctr-maxrng/2
        zmin = zctr-maxrng/2
        xmax = xctr+maxrng/2
        ymax = yctr+maxrng/2
        zmax = zctr+maxrng/2
        # x = [xmin, xmin, xmin, xmin, xmax, xmax, xmax, xmax]
        # y = [ymin, ymin, ymax, ymax, ymin, ymin, ymax, ymax]
        # z = [zmin, zmax, zmin, zmax, zmin, zmax, zmin, zmax]
        # ax.scatter(x, y, z)
        ax.set_xlim3d(xmin, xmax)
        ax.set_ylim3d(ymin, ymax)
        ax.set_zlim3d(zmin, zmax)

def set_axes_radius(ax, origin, radius):
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

def set_axes_equal(ax):
    from numpy import array, mean, max, abs
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    limits = array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = mean(limits, axis=1)
    radius = 0.5 * max(abs(limits[:, 1] - limits[:, 0]))
    set_axes_radius(ax, origin, radius)
