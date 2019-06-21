from pygeom.geom3d import Vector, jhat
from numpy.matlib import zeros, empty, matrix
from math import pi, radians, cos, sin
from matplotlib.pyplot import figure
from .latticesystem import LatticeSystem

class LatticeResult(object):
    name = None
    sys = None
    rho = None
    alpha = None
    beta = None
    speed = None
    _vfs = None
    _qfs = None
    _udc = None
    _ulc = None
    _uyc = None
    _gam = None
    _gmat = None
    _avv = None
    _afv = None
    _amv = None
    _phi = None
    _pmat = None
    _nfvel = None
    _nffrc = None
    _nfmom = None
    _nffrctot = None
    _nfmomtot = None
    _trlft = None
    _trdrg = None
    _trfrc = None
    _trwsh = None
    _strpy = None
    _Cx = None
    _Cy = None
    _Cz = None
    _Cl = None
    _Cm = None
    _Cn = None
    _CL = None
    _CDi = None
    _CY = None
    _CL_ff = None
    _CDi_ff = None
    _CY_ff = None
    _e = None
    def __init__(self, name: str, sys: LatticeSystem):
        self.name = name
        self.sys = sys
    def reset(self):
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None
    def set_conditions(self, alpha: float=0.0, beta: float=0.0, speed: float=1.0, rho: float=1.0):
        self.alpha = alpha
        self.beta = beta
        self.speed = speed
        self.rho = rho
    def set_freestream_velocity(self, vfs: Vector):
        self.vfs = vfs
        self.udc = self.vfs.to_unit()
        self.ulc = (self.udc**jhat).to_unit()
        self.uyc = (self.udc**self.ulc).to_unit()
        num = len(self.sys.pnls)
        self.avv = empty((num, 1), dtype=Vector)
        self.afv = empty((num, 1), dtype=Vector)
        self.amv = empty((num, 1), dtype=Vector)
        for i, pnl in enumerate(self.sys.pnls):
            self.avv[i, 0] = self.vfs
            self.afv[i, 0] = self.avv[i, 0]**pnl.leni
            self.amv[i, 0] = (pnl.pnti-self.sys.rref)**self.afv[i, 0]
    def set_gam(self, gam: list):
        self._gam = gam
        self._gmat = matrix(self._gam, dtype=float).transpose()
    def set_phi(self, phi: list):
        if len(phi) != len(self.sys.strps):
            raise Exception('The length of phi must equal the number of strips.')
        self._phi = phi
        self._pmat = matrix(self._phi, dtype=float).transpose()
    def set_lift_distribution(self, l: list, rho: float, speed: float):
        if len(l) != len(self.sys.strps):
            raise Exception('The length of l must equal the number of strips.')
        self._phi = [li/rho/speed for li in l]
        self._pmat = matrix(self._phi, dtype=float).transpose()
    @property
    def gmat(self):
        if self._gmat is None:
            gamx = self.sys.gam[:, 0]
            gamy = self.sys.gam[:, 1]
            gamz = self.sys.gam[:, 2]
            self._gmat = gamx*self.vfs.x+gamy*self.vfs.y+gamz*self.vfs.z
        return self._gmat
    @property
    def gam(self):
        if self._gam is None:
            self._gam = [self.gmat[i, 0] for i in range(self.gmat.shape[0])]
        return self._gam
    @property
    def strpy(self):
        if self._strpy is None:
            self._strpy = [strp.pnti.y for strp in self.sys.strps]
        return self._strpy
    @property
    def pmat(self):
        if self._pmat is None:
            num = len(self.sys.strps)
            self._pmat = zeros((num, 1))
            for i, strp in enumerate(self.sys.strps):
                for pnl in strp.pnls:
                    self._pmat[i, 0] += self.gam[pnl.lpid]
        return self._pmat
    @property
    def phi(self):
        if self._phi is None:
            self._phi = [self.pmat[i, 0] for i in range(self.pmat.shape[0])]
        return self._phi
    @property
    def vfs(self):
        if self._vfs is None:
            if self.alpha is None:
                self.alpha = 0.0
            if self.beta is None:
                self.beta = 0.0
            if self.speed is None:
                self.speed = 1.0
            alrad = radians(self.alpha)
            btrad = radians(self.beta)
            cosal = cos(alrad)
            sinal = sin(alrad)
            cosbt = cos(btrad)
            sinbt = sin(btrad)
            self._vfs = Vector(cosal*cosbt, -sinbt, sinal*cosbt)*self.speed
        return self._vfs
    @property
    def qfs(self):
        if self._qfs is None:
            self._qfs = self.rho*self.speed**2/2
        return self._qfs
    @property
    def udc(self):
        if self._udc is None:
            self._udc = self.vfs.to_unit()
        return self._udc
    @property
    def ulc(self):
        if self._ulc is None:
            self._ulc = (self.udc**jhat).to_unit()
        return self._ulc
    @property
    def uyc(self):
        if self._uyc is None:
            self._uyc = (self.ulc**self.udc).to_unit()
        return self._uyc
    @property
    def avv(self):
        if self._avv is None:
            num = len(self.sys.pnls)
            self._avv = empty((num, 1), dtype=Vector)
            for i in range(num):
                self._avv[i, 0] = self.vfs
        return self._avv
    @property
    def afv(self):
        if self._afv is None:
            num = len(self.sys.pnls)
            self._afv = empty((num, 1), dtype=Vector)
            for i, pnl in enumerate(self.sys.pnls):
                self._afv[i, 0] = self.avv[i, 0]**pnl.leni
        return self._afv
    @property
    def amv(self):
        if self._amv is None:
            num = len(self.sys.pnls)
            self._amv = empty((num, 1), dtype=Vector)
            for i, pnl in enumerate(self.sys.pnls):
                self._amv[i, 0] = (pnl.pnti-self.sys.rref)**self.afv[i, 0]
        return self._amv
    @property
    def nfvel(self):
        if self._nfvel is None:
            tmp = self.sys.avg*self.gmat+self.avv
            self._nfvel = [tmp[i, 0] for i in range(tmp.shape[0])]
        return self._nfvel
    @property
    def nffrc(self):
        if self._nffrc is None:
            tmp = self.sys.afg*self.gmat+self.afv
            self._nffrc = [gam*tmp[i, 0] for i, gam in enumerate(self.gam)]
        return self._nffrc
    @property
    def nfmom(self):
        if self._nfmom is None:
            tmp = self.sys.amg*self.gmat+self.amv
            self._nfmom = [gam*tmp[i, 0] for i, gam in enumerate(self.gam)]
        return self._nfmom
    @property
    def trlft(self):
        if self._trlft is None:
            rhoV = self.rho*self.speed
            self._trlft = [phi*self.sys.blg[i, 0]*rhoV for i, phi in enumerate(self.phi)]
        return self._trlft
    @property
    def trdrg(self):
        if self._trdrg is None:
            rho = self.rho
            tmp = self.sys.bdg*self.pmat
            self._trdrg = [phi*tmp[i, 0]*rho for i, phi in enumerate(self.phi)]
        return self._trdrg
    @property
    def trfrc(self):
        if self._trfrc is None:
            rhoV = self.rho*self.speed
            self._trfrc = [phi*self.sys.byg[i, 0]*rhoV for i, phi in enumerate(self.phi)]
        return self._trfrc
    @property
    def trwsh(self):
        if self._trwsh is None:
            tmp = self.sys.bvg*self.pmat
            self._trwsh = [tmp[i, 0] for i in range(tmp.shape[0])]
        return self._trwsh
    @property
    def nffrctot(self):
        if self._nffrctot is None:
            self._nffrctot = sum(self.nffrc)
        return self._nffrctot
    @property
    def nfmomtot(self):
        if self._nfmomtot is None:
            self._nfmomtot = sum(self.nfmom)
        return self._nfmomtot
    @property
    def Cx(self):
        if self._Cx is None:
            self._Cx = -self.nffrctot.x/self.qfs/self.sys.sref
        return self._Cx
    @property
    def Cy(self):
        if self._Cy is None:
            self._Cy = self.nffrctot.y/self.qfs/self.sys.sref
        return self._Cy
    @property
    def Cz(self):
        if self._Cz is None:
            self._Cz = -self.nffrctot.z/self.qfs/self.sys.sref
        return self._Cz
    @property
    def CL(self):
        if self._CL is None:
            self._CL = self.ulc*self.nffrctot/self.qfs/self.sys.sref
        return self._CL
    @property
    def CDi(self):
        if self._CDi is None:
            self._CDi = self.udc*self.nffrctot/self.qfs/self.sys.sref
        return self._CDi
    @property
    def CY(self):
        if self._CY is None:
            self._CY = self.uyc*self.nffrctot/self.qfs/self.sys.sref
        return self._CY
    @property
    def Cl(self):
        if self._Cl is None:
            self._Cl = -self.nfmomtot.x/self.qfs/self.sys.sref/self.sys.bref
        return self._Cl
    @property
    def Cm(self):
        if self._Cm is None:
            self._Cm = self.nfmomtot.y/self.qfs/self.sys.sref/self.sys.cref
        return self._Cm
    @property
    def Cn(self):
        if self._Cn is None:
            self._Cn = -self.nfmomtot.z/self.qfs/self.sys.sref/self.sys.bref
        return self._Cn
    @property
    def CL_ff(self):
        if self._CL_ff is None:
            self._CL_ff = sum(self.trlft)/self.qfs/self.sys.sref
        return self._CL_ff
    @property
    def CDi_ff(self):
        if self._CDi_ff is None:
            self._CDi_ff = sum(self.trdrg)/self.qfs/self.sys.sref
        return self._CDi_ff
    @property
    def CY_ff(self):
        if self._CY_ff is None:
            self._CY_ff = sum(self.trfrc)/self.qfs/self.sys.sref
        return self._CY_ff
    @property
    def e(self):
        if self._e is None:
            self._e = self.CL_ff**2/pi/self.sys.ar/self.CDi_ff
        return self._e
    def plot_phi_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.strpy, self.phi, label=self.name)
        ax.legend()
        return ax
    def plot_trefftz_lift_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        trlft = [self.trlft[i]/strp.dyt for i, strp in enumerate(self.sys.strps)]
        ax.plot(self.strpy, trlft, label=self.name)
        ax.legend()
        return ax
    def plot_trefftz_drag_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        trdrg = [self.trdrg[i]/strp.dst for i, strp in enumerate(self.sys.strps)]
        ax.plot(self.strpy, trdrg, label=self.name)
        ax.legend()
        return ax
    def plot_trefftz_wash_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.strpy, self.trwsh, label=self.name)
        ax.legend()
        return ax
    def print_near_field_total_loads(self):
        print(f'Total Force in X = {self.nffrctot.x}')
        print(f'Total Force in Y = {self.nffrctot.y}')
        print(f'Total Force in Z = {self.nffrctot.z}')
        print(f'Total Moment in X = {self.nfmomtot.x}')
        print(f'Total Moment in Y = {self.nfmomtot.y}')
        print(f'Total Moment in Z = {self.nfmomtot.z}')
    def print_aerodynamic_coefficients(self):
        from . import cfrm, dfrm
        if self._gam is not None:
            print('Cx = '+cfrm.format(self.Cx))
            print('Cy = '+cfrm.format(self.Cy))
            print('Cz = '+cfrm.format(self.Cz))
            print('Cl = '+cfrm.format(self.Cl))
            print('Cm = '+cfrm.format(self.Cm))
            print('Cn = '+cfrm.format(self.Cn))
            print('CL = '+cfrm.format(self.CL))
            print('CDi = '+dfrm.format(self.CDi))
            print('CY = '+cfrm.format(self.CY))
        if self._phi is not None:
            print('CL_ff = '+cfrm.format(self.CL_ff))
            print('CDi_ff = '+dfrm.format(self.CDi_ff))
            print('CY_ff = '+cfrm.format(self.CY_ff))
            print('e = '+cfrm.format(self.e))
