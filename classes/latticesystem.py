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
    def optimum_lift_distribution(self, Lspec: float, lspec: float=None):
        # from IPython.display import display
        # from pyhtml import HTMLMatrix
        nump = len(self.strps)
        numc = 1
        if lspec is not None:
            numc += 2
        amat = zeros((nump+numc, nump+numc))
        bmat = zeros((nump+numc, 1))
        amat[0:nump, 0:nump] = self.bdg + self.bdg.transpose()
        amat[0:nump, nump] = self.blg
        amat[nump, 0:nump] = self.blg.transpose()
        bmat[nump, 0] = Lspec
        bmg1 = self.bmg.copy()
        bmg2 = self.bmg.copy()
        hlfp = int(nump/2)
        for i in range(hlfp):
            bmg1[i, 0] = 0.0
        for i in range(hlfp, nump):
            bmg2[i, 0] = 0.0
        if lspec is not None:
            amat[0:nump, nump+1] = bmg1
            amat[nump+1, 0:nump] = bmg1.transpose()
            bmat[nump+1, 0] = lspec
            amat[0:nump, nump+2] = bmg2
            amat[nump+2, 0:nump] = bmg2.transpose()
            bmat[nump+2, 0] = -1.0*lspec
        xmat = solve(amat, bmat)
        phi = [xmat[i, 0] for i in range(nump)]
        lam = [xmat[nump+j, 0] for j in range(numc)]
        Di = (xmat[0:nump, 0].transpose()*self.bdg*xmat[0:nump, 0])[0, 0]
        L = (xmat[0:nump, 0].transpose()*self.blg)[0, 0]
        l = (xmat[0:nump, 0].transpose()*bmg1)[0, 0]
        # amat = HTMLMatrix(amat)
        # display(amat)
        return phi, lam, Di, L, l
    @property
    def ar(self):
        if self._ar is None:
            self._ar = self.bref**2/self.sref
        return self._ar
    def set_strip_alpha(self, alpha: list):
        for i, strp in enumerate(self.strps):
            strp.alpha = alpha[i]
        self._aic = None
        self._afs = None
        self._gam = None
    def __repr__(self):
        return '<LatticeSystem {:s}>'.format(self.name)

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
