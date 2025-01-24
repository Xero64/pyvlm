from typing import TYPE_CHECKING, Any

from matplotlib.pyplot import figure
from numpy import cos, pi, radians, sin, tan, zeros
from py2md.classes import MDReport
from pygeom.geom3d import Coordinate, Vector
from ..tools.mass import Mass

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from mpl_toolkits.mplot3d import Axes3D
    from numpy.typing import NDArray

    from .latticesystem import LatticeSystem as System


class LatticeResult():
    name: str = None
    sys: 'System' = None
    rho: float = None
    mach: float = None
    speed: float = None
    alpha: float = None
    beta: float = None
    pbo2V: float = None
    qco2V: float = None
    rbo2V: float = None
    ctrls: dict[str, float] = None
    rcg: Vector = None
    mass: Mass = None
    _acs: Coordinate = None
    _scs: Coordinate = None
    _dacsa: dict[str, Vector] = None
    _dacsb: dict[str, Vector] = None
    _dscsa: dict[str, Vector] = None
    _dscsb: dict[str, Vector] = None
    _vfs: Vector = None
    _dvfsa: Vector = None
    _dvfsb: Vector = None
    _pqr: Vector = None
    _ofs: Vector = None
    _qfs: float = None
    _arm: Vector = None
    _ungam: Vector = None
    _gamma: 'NDArray' = None
    _avg: Vector = None
    _avv: Vector = None
    _afg: Vector = None
    _afv: Vector = None
    _phi: 'NDArray' = None
    _bvv: Vector = None
    _brm: Vector = None
    _nfres: 'GammaResult' = None
    _trres: 'PhiResult' = None
    _pdres: 'ParasiticDragResult' = None
    _stripres: 'StripResult' = None
    _stres: 'StabilityResult' = None
    _ctgamp: dict[str, 'NDArray'] = None
    _ctresp: dict[str, 'GammaResult'] = None
    _ctgamn: dict[str, 'NDArray'] = None
    _ctresn: dict[str, 'GammaResult'] = None

    def __init__(self, name: str, sys: 'System') -> None:
        self.name = name
        self.sys = sys
        self.initialise()

    def initialise(self) -> None:
        self.rho = 1.0
        self.mach = 0.0
        self.speed = 1.0
        self.alpha = 0.0
        self.beta = 0.0
        self.pbo2V = 0.0
        self.qco2V = 0.0
        self.rbo2V = 0.0
        self.ctrls = {}
        for control in self.sys.ctrls:
            self.ctrls[control] = 0.0
        self.rcg = self.sys.rref

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr.startswith('_'):
                setattr(self, attr, None)

    def set_density(self, rho: float) -> None:
        self.rho = rho
        self.reset()

    def set_state(self, mach: float | None = None, speed: float | None = None,
                  alpha: float | None = None, beta: float | None = None,
                  pbo2V: float | None = None, qco2V: float | None = None, rbo2V: float | None = None) -> None:
        if mach is not None:
            self.mach = mach
        if speed is not None:
            self.speed = speed
        if alpha is not None:
            self.alpha = alpha
        if beta is not None:
            self.beta = beta
        if pbo2V is not None:
            self.pbo2V = pbo2V
        if qco2V is not None:
            self.qco2V = qco2V
        if rbo2V is not None:
            self.rbo2V = rbo2V
        self.reset()

    def set_controls(self, **kwargs: dict[str, float]) -> None:
        for control in kwargs:
            self.ctrls[control] = kwargs[control]
        self.reset()

    def set_cg(self, rcg: Vector) -> None:
        self.rcg = rcg
        self.reset()

    def calc_coordinate_systems(self) -> None:
        pnt = self.rcg
        cosal, sinal = trig_angle(self.alpha)
        cosbt, sinbt = trig_angle(self.beta)
        self._vfs = Vector(cosal*cosbt, -sinbt, sinal*cosbt)*self.speed
        self._dvfsa = Vector(-sinal*cosbt, 0.0, cosal*cosbt)*self.speed
        self._dvfsb = Vector(-cosal*sinbt, -cosbt, -sinal*sinbt)*self.speed
        # Aerodynamic Coordinate System
        # acs_dirx = Vector(cosal*cosbt, -sinbt, sinal*cosbt)
        # acs_diry = Vector(cosal*sinbt, cosbt, sinal*sinbt)
        acs_dirx = Vector(cosal, 0.0, sinal)
        acs_diry = Vector(0.0, 1.0, 0.0)
        self._acs = Coordinate(pnt, acs_dirx, acs_diry)
        # Stability Coordinate System
        scs_dirx = Vector(-cosal, 0.0, -sinal)
        scs_diry = Vector(0.0, 1.0, 0.0)
        self._scs = Coordinate(pnt, scs_dirx, scs_diry)
        # Derivative of Aerodynamic Coordinate System wrt alpha
        # self._dacsa = {'x': Vector(-sinal*cosbt, 0.0, cosal*cosbt),
        #                'y': Vector(-sinal*sinbt, 0.0, cosal*sinbt),
        #                'z': Vector(-cosal, 0.0, -sinal)}
        self._dacsa = {'x': Vector(-sinal, 0.0, cosal),
                       'y': Vector(0.0, 0.0, 0.0),
                       'z': Vector(-cosal, 0.0, -sinal)}
        # Derivative of Aerodynamic Coordinate System wrt beta
        # self._dacsb = {'x': Vector(-cosal*sinbt, -cosbt, -sinal*sinbt),
        #                'y': Vector(cosal*cosbt, -sinbt, sinal*cosbt),
        #                'z': Vector(0.0, 0.0, 0.0)}
        self._dacsb = {'x': Vector.zeros(),
                       'y': Vector.zeros(),
                       'z': Vector.zeros()}
        # Derivative of Stability Coordinate System wrt alpha
        self._dscsa = {'x': Vector(sinal, 0.0, -cosal),
                       'y': Vector(0.0, 0.0, 0.0),
                       'z': Vector(cosal, 0.0, sinal)}
        # Derivative of Stability Coordinate System wrt beta
        self._dscsb = {'x': Vector.zeros(),
                       'y': Vector.zeros(),
                       'z': Vector.zeros()}

    @property
    def acs(self) -> Coordinate:
        if self._acs is None:
            self.calc_coordinate_systems()
        return self._acs

    @property
    def scs(self) -> Coordinate:
        if self._scs is None:
            self.calc_coordinate_systems()
        return self._scs

    @property
    def dacsa(self) -> dict[str, Vector]:
        if self._dacsa is None:
            self.calc_coordinate_systems()
        return self._dacsa

    @property
    def dacsb(self) -> dict[str, Vector]:
        if self._dacsb is None:
            self.calc_coordinate_systems()
        return self._dacsb

    @property
    def dscsa(self) -> dict[str, Vector]:
        if self._dscsa is None:
            self.calc_coordinate_systems()
        return self._dscsa

    @property
    def dscsb(self) -> dict[str, Vector]:
        if self._dscsb is None:
            self.calc_coordinate_systems()
        return self._dscsb

    @property
    def vfs(self) -> Vector:
        if self._vfs is None:
            self.calc_coordinate_systems()
        return self._vfs

    @property
    def dvfsa(self) -> Vector:
        if self._dvfsa is None:
            self.calc_coordinate_systems()
        return self._dvfsa

    @property
    def dvfsb(self) -> Vector:
        if self._dvfsb is None:
            self.calc_coordinate_systems()
        return self._dvfsb

    @property
    def pqr(self) -> Vector:
        if self._pqr is None:
            p = self.pbo2V*2*self.speed/self.sys.bref
            q = self.qco2V*2*self.speed/self.sys.cref
            r = self.rbo2V*2*self.speed/self.sys.bref
            self._pqr = Vector(p, q, r)
        return self._pqr

    @property
    def ofs(self) -> Vector:
        if self._ofs is None:
            self._ofs = self.scs.vector_to_global(self.pqr)
        return self._ofs

    @property
    def qfs(self) -> float:
        if self._qfs is None:
            self._qfs = self.rho*self.speed**2/2
        return self._qfs

    @property
    def arm(self) -> Vector:
        if self._arm is None:
            num = len(self.sys.pnls)
            self._arm = Vector.zeros(num)
            for pnl in self.sys.pnls:
                i = pnl.lpid
                self._arm[i] = pnl.pnti - self.rcg
        return self._arm

    @property
    def ungam(self) -> Vector:
        if self._ungam is None:
            self._ungam = self.sys.ungam(self.mach)
        return self._ungam

    @property
    def gamma(self):
        if self._gamma is None:
            self._gamma = self.ungam[:, 0].dot(self.vfs)
            self._gamma += self.ungam[:, 1].dot(self.ofs)
            for control in self.ctrls:
                if control in self.sys.ctrls:
                    ctrl = self.ctrls[control]
                    ctrlrad = radians(ctrl)
                    index = self.sys.ctrls[control]
                    if ctrl >= 0.0:
                        indv = index[0]
                        indo = index[1]
                    else:
                        indv = index[2]
                        indo = index[3]
                    self._gamma += ctrlrad*(self.ungam[:, indv].dot(self.vfs))
                    self._gamma += ctrlrad*(self.ungam[:, indo].dot(self.ofs))
        return self._gamma

    def gctrlp_single(self, control: str) -> 'NDArray':
        indv = self.sys.ctrls[control][0]
        indo = self.sys.ctrls[control][1]
        return self.ungam[:, indv].dot(self.vfs) + self.ungam[:, indo].dot(self.ofs)

    def gctrlp(self) -> dict[str, 'NDArray']:
        gmats = {}
        for control in self.sys.ctrls:
            gmats[control] = self.gctrlp_single(control)
        return gmats

    def gctrln_single(self, control: str) -> 'NDArray':
        indv = self.sys.ctrls[control][2]
        indo = self.sys.ctrls[control][3]
        return self.ungam[:, indv].dot(self.vfs) + self.ungam[:, indo].dot(self.ofs)

    def gctrln(self) -> dict[str, 'NDArray']:
        gmats = {}
        for control in self.sys.ctrls:
            gmats[control] = self.gctrln_single(control)
        return gmats

    @property
    def phi(self) -> 'NDArray':
        if self._phi is None:
            num = len(self.sys.strps)
            self._phi = zeros(num)
            for strp in self.sys.strps:
                for pnl in strp.pnls:
                    self._phi[strp.lsid] += self.gamma[pnl.lpid]
        return self._phi

    @property
    def avg(self) -> Vector:
        if self._avg is None:
            self._avg = self.sys.avg(self.mach)
        return self._avg

    @property
    def avv(self) -> Vector:
        if self._avv is None:
            num = len(self.sys.pnls)
            self._avv = Vector.zeros(num)
            for pnl in self.sys.pnls:
                i = pnl.lpid
                self._avv[i] = self.vfs - self.ofs.cross(self.arm[i])
        return self._avv

    @property
    def bvv(self) -> Vector:
        if self._bvv is None:
            num = len(self.sys.strps)
            self._bvv = Vector.zeros(num)
            for strp in self.sys.strps:
                if not strp.noload:
                    i = strp.lsid
                    self._bvv[i] = self.vfs - self.ofs.cross(self.brm[i])
        return self._bvv

    @property
    def brm(self) -> Vector:
        if self._brm is None:
            num = len(self.sys.strps)
            self._brm = Vector.zeros(num)
            for strp in self.sys.strps:
                i = strp.lsid
                self._brm[i] = strp.pntq - self.rcg
        return self._brm

    @property
    def afg(self) -> Vector:
        if self._afg is None:
            self._afg = self.sys.afg(self.mach)
        return self._afg

    @property
    def afv(self) -> Vector:
        if self._afv is None:
            num = len(self.sys.pnls)
            self._afv = Vector.zeros(num)
            for pnl in self.sys.pnls:
                if not pnl.noload:
                    i = pnl.lpid
                    self._afv[i] = self.avv[i].cross(pnl.leni)
        return self._afv

    @property
    def nfres(self) -> 'GammaResult':
        if self._nfres is None:
            self._nfres = GammaResult(self, self.gamma)
        return self._nfres

    @property
    def trres(self) -> 'PhiResult':
        if self._trres is None:
            self._trres = PhiResult(self, self.phi)
        return self._trres

    @property
    def pdres(self) -> 'ParasiticDragResult':
        if self._pdres is None:
            self._pdres = ParasiticDragResult(self)
        return self._pdres

    @property
    def stripres(self) -> 'StripResult':
        if self._stripres is None:
            self._stripres = StripResult(self.nfres)
        return self._stripres

    @property
    def stres(self) -> 'StabilityResult':
        if self._stres is None:
            self._stres = StabilityResult(self)
        return self._stres

    @property
    def ctgamp(self):
        if self._ctgamp is None:
            self._ctgamp = self.gctrlp()
        return self._ctgamp

    @property
    def ctresp(self):
        if self._ctresp is None:
            self._ctresp = {}
            for control in self.ctgamp:
                self._ctresp[control] = GammaResult(self, self.ctgamp[control])
        return self._ctresp

    @property
    def ctgamn(self):
        if self._ctgamn is None:
            self._ctgamn = self.gctrln()
        return self._ctgamn

    @property
    def ctresn(self):
        if self._ctresn is None:
            self._ctresn = {}
            for control in self.ctgamn:
                self._ctresn[control] = GammaResult(self, self.ctgamn[control])
        return self._ctresn

    def plot_panel_near_field_velocities(self, ax: 'Axes' = None,
                                         component: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        py = [pnl.pnti.y for pnl in self.sys.pnls]
        if component is None or component == 'x':
            ax.plot(py, self.nfres.nfvel.x, label=self.name+' Velocity X')
        if component is None or component == 'y':
            ax.plot(py, self.nfres.nfvel.y, label=self.name+' Velocity Y')
        if component is None or component == 'z':
            ax.plot(py, self.nfres.nfvel.z, label=self.name+' Velocity Z')
        ax.legend()
        return ax

    def plot_phi_distribution(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        for srfc in self.sys.srfcs:
            ax.plot(srfc.strpy, self.phi[srfc.strpi], label=self.name+' for '+srfc.name)
        ax.legend()
        return ax

    def plot_strip_lift_force_distribution(self, ax: 'Axes' = None, axis: str='b',
                                           surfaces: list | None = None, normalise: bool=False,
                                           label: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        if surfaces is None:
            srfcs = self.sys.srfcs
        else:
            srfcs = []
            for srfc in self.sys.srfcs:
                if srfc.name in surfaces:
                    srfcs.append(srfc)
        onesrfc = len(srfcs) == 1
        for srfc in srfcs:
            if normalise:
                l = [self.stripres.lift[strp.lsid]/strp.area/self.qfs for strp in srfc.strps]
            else:
                l = [self.stripres.lift[strp.lsid]/strp.dst for strp in srfc.strps]
            if label is None:
                thislabel = self.name+' for '+srfc.name
            else:
                if not onesrfc:
                    thislabel = label+' for '+srfc.name
                else:
                    thislabel = label
            if axis == 'b':
                b = srfc.strpb
                ax.plot(b, l, label=thislabel)
            elif axis == 'y':
                y = srfc.strpy
                if max(y) > min(y):
                    ax.plot(y, l, label=thislabel)
            elif axis == 'z':
                z = srfc.strpz
                if max(z) > min(z):
                    ax.plot(l, z, label=thislabel)
        ax.legend()
        return ax

    def plot_strip_side_force_distribution(self, ax: 'Axes' = None, axis: str='b',
                                           surfaces: list | None = None, normalise: bool=False,
                                           label: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        if surfaces is None:
            srfcs = self.sys.srfcs
        else:
            srfcs = []
            for srfc in self.sys.srfcs:
                if srfc.name in surfaces:
                    srfcs.append(srfc)
        onesrfc = len(srfcs) == 1
        for srfc in srfcs:
            if normalise:
                f = [self.stripres.side[strp.lsid]/strp.area/self.qfs for strp in srfc.strps]
            else:
                f = [self.stripres.side[strp.lsid]/strp.dst for strp in srfc.strps]
            if label is None:
                thislabel = self.name+' for '+srfc.name
            else:
                if not onesrfc:
                    thislabel = label+' for '+srfc.name
                else:
                    thislabel = label
            if axis == 'b':
                b = srfc.strpb
                ax.plot(b, f, label=thislabel)
            elif axis == 'y':
                y = srfc.strpy
                if max(y) > min(y):
                    ax.plot(y, f, label=thislabel)
            elif axis == 'z':
                z = srfc.strpz
                if max(z) > min(z):
                    ax.plot(f, z, label=thislabel)
        ax.legend()
        return ax

    def plot_strip_drag_force_distribution(self, ax: 'Axes' = None, axis: str='b',
                                           surfaces: list | None = None, normalise: bool=False,
                                           label: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        if surfaces is None:
            srfcs = self.sys.srfcs
        else:
            srfcs = []
            for srfc in self.sys.srfcs:
                if srfc.name in surfaces:
                    srfcs.append(srfc)
        onesrfc = len(srfcs) == 1
        for srfc in srfcs:
            if normalise:
                d = [self.stripres.drag[strp.lsid]/strp.area/self.qfs for strp in srfc.strps]
            else:
                d = [self.stripres.drag[strp.lsid]/strp.dst for strp in srfc.strps]
            if label is None:
                thislabel = self.name+' for '+srfc.name
            else:
                if not onesrfc:
                    thislabel = label+' for '+srfc.name
                else:
                    thislabel = label
            if axis == 'b':
                b = srfc.strpb
                ax.plot(b, d, label=thislabel)
            elif axis == 'y':
                y = srfc.strpy
                if max(y) > min(y):
                    ax.plot(y, d, label=thislabel)
            elif axis == 'z':
                z = srfc.strpz
                if max(z) > min(z):
                    ax.plot(d, z, label=thislabel)
        ax.legend()
        return ax

    def plot_trefftz_lift_force_distribution(self, ax: 'Axes' = None, axis: str='b',
                                             surfaces: list | None = None, normalise: bool=False,
                                             label: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        if surfaces is None:
            srfcs = self.sys.srfcs
        else:
            srfcs = []
            for srfc in self.sys.srfcs:
                if srfc.name in surfaces:
                    srfcs.append(srfc)
        onesrfc = len(srfcs) == 1
        for srfc in srfcs:
            if normalise:
                l = [self.trres.trfrc.z[strp.lsid]/strp.area/self.qfs for strp in srfc.strps]
            else:
                l = [self.trres.trfrc.z[strp.lsid]/strp.dst for strp in srfc.strps]
            if label is None:
                thislabel = self.name+' for '+srfc.name
            else:
                if not onesrfc:
                    thislabel = label+' for '+srfc.name
                else:
                    thislabel = label
            if axis == 'b':
                b = srfc.strpb
                ax.plot(b, l, label=thislabel)
            elif axis == 'y':
                y = srfc.strpy
                if max(y) > min(y):
                    ax.plot(y, l, label=thislabel)
            elif axis == 'z':
                z = srfc.strpz
                if max(z) > min(z):
                    ax.plot(l, z, label=thislabel)
        ax.legend()
        return ax

    def plot_trefftz_side_force_distribution(self, ax: 'Axes' = None, axis: str='b',
                                             surfaces: list | None = None, normalise: bool=False,
                                             label: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        if surfaces is None:
            srfcs = self.sys.srfcs
        else:
            srfcs = []
            for srfc in self.sys.srfcs:
                if srfc.name in surfaces:
                    srfcs.append(srfc)
        onesrfc = len(srfcs) == 1
        for srfc in srfcs:
            if normalise:
                f = [self.trres.trfrc.y[strp.lsid]/strp.area/self.qfs for strp in srfc.strps]
            else:
                f = [self.trres.trfrc.y[strp.lsid]/strp.dst for strp in srfc.strps]
            if label is None:
                thislabel = self.name+' for '+srfc.name
            else:
                if not onesrfc:
                    thislabel = label+' for '+srfc.name
                else:
                    thislabel = label
            if axis == 'b':
                b = srfc.strpb
                ax.plot(b, f, label=thislabel)
            elif axis == 'y':
                y = srfc.strpy
                if max(y) > min(y):
                    ax.plot(y, f, label=thislabel)
            elif axis == 'z':
                z = srfc.strpz
                if max(z) > min(z):
                    ax.plot(f, z, label=thislabel)
        ax.legend()
        return ax

    def plot_trefftz_drag_force_distribution(self, ax: 'Axes' = None, axis: str='b',
                                             surfaces: list | None = None, normalise: bool=False,
                                             label: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        if surfaces is None:
            srfcs = self.sys.srfcs
        else:
            srfcs = []
            for srfc in self.sys.srfcs:
                if srfc.name in surfaces:
                    srfcs.append(srfc)
        onesrfc = len(srfcs) == 1
        for srfc in srfcs:
            if normalise:
                d = [self.trres.trfrc.x[strp.lsid]/strp.area/self.qfs for strp in srfc.strps]
            else:
                d = [self.trres.trfrc.x[strp.lsid]/strp.dst for strp in srfc.strps]
            if label is None:
                thislabel = self.name+' for '+srfc.name
            else:
                if not onesrfc:
                    thislabel = label+' for '+srfc.name
                else:
                    thislabel = label
            if axis == 'b':
                b = srfc.strpb
                ax.plot(b, d, label=thislabel)
            elif axis == 'y':
                y = srfc.strpy
                if max(y) > min(y):
                    ax.plot(y, d, label=thislabel)
            elif axis == 'z':
                z = srfc.strpz
                if max(z) > min(z):
                    ax.plot(d, z, label=thislabel)
        ax.legend()
        return ax

    def plot_trefftz_wash_distribution(self, ax: 'Axes' = None, axis: str='b',
                                       surfaces: list | None = None, normalise: bool=False,
                                       label: str | None = None) -> 'Axes':
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        if surfaces is None:
            srfcs = self.sys.srfcs
        else:
            srfcs = []
            for srfc in self.sys.srfcs:
                if srfc.name in surfaces:
                    srfcs.append(srfc)
        onesrfc = len(srfcs) == 1
        for srfc in srfcs:
            if normalise:
                w = [self.trres.trwsh[strp.lsid]/self.speed for strp in srfc.strps]
            else:
                w = [self.trres.trwsh[strp.lsid] for strp in srfc.strps]
            if label is None:
                thislabel = self.name+' for '+srfc.name
            else:
                if not onesrfc:
                    thislabel = label+' for '+srfc.name
                else:
                    thislabel = label
            if axis == 'b':
                b = srfc.strpb
                ax.plot(b, w, label=thislabel)
            elif axis == 'y':
                y = srfc.strpy
                if max(y) > min(y):
                    ax.plot(y, w, label=thislabel)
            elif axis == 'z':
                z = srfc.strpz
                if max(z) > min(z):
                    ax.plot(w, z, label=thislabel)
        ax.legend()
        return ax

    def plot_surfaces(self) -> 'Axes3D':
        ax = None
        for srfc in self.sys.srfcs:
            ax = srfc.plot_surface(ax=ax)
        return

    def to_result(self, name: str='') -> 'LatticeResult':
        if name == '':
            name = self.name
        res = LatticeResult(name, self.sys)
        res.set_density(rho=self.rho)
        res.set_state(speed=self.speed, alpha=self.alpha, beta=self.beta,
                      pbo2V=self.pbo2V, qco2V=self.qco2V, rbo2V=self.rbo2V)
        res.set_controls(**self.ctrls)
        res.set_cg(self.rcg)
        return res

    @property
    def surface_loads(self) -> MDReport:
        report = MDReport()
        report.add_heading('Surface Loads', 2)
        table = report.add_table()
        table.add_column('xref', '.3f', data=[self.rcg.x])
        table.add_column('yref', '.3f', data=[self.rcg.y])
        table.add_column('zref', '.3f', data=[self.rcg.z])
        table1 = report.add_table()
        table1.add_column('Name', 's')
        table1.add_column('Fx', '.3f')
        table1.add_column('Fy', '.3f')
        table1.add_column('Fz', '.3f')
        table1.add_column('Mx', '.3f')
        table1.add_column('My', '.3f')
        table1.add_column('Mz', '.3f')
        table2 = report.add_table()
        table2.add_column('Name', 's')
        table2.add_column('Area', '.3f')
        table2.add_column('Di', '.3f')
        table2.add_column('Y', '.3f')
        table2.add_column('L', '.3f')
        table2.add_column('CDi', '.7f')
        table2.add_column('CY', '.5f')
        table2.add_column('CL', '.5f')
        Ditot = 0.0
        Ytot = 0.0
        Ltot = 0.0
        for srfc in self.sys.srfcs:
            area = srfc.area
            ind = srfc.pnli
            frc = self.nfres.nffrc[ind].sum()
            mom = self.nfres.nfmom[ind].sum()
            table1.add_row([srfc.name, frc.x, frc.y, frc.z, mom.x, mom.y, mom.z])
            if area > 0.0:
                Di = frc.dot(self.acs.dirx)
                Y = frc.dot(self.acs.diry)
                L = frc.dot(self.acs.dirz)
                CDi = Di/self.qfs/area
                CY = Y/self.qfs/area
                CL = L/self.qfs/area
                table2.add_row([srfc.name, area, Di, Y, L, CDi, CY, CL])
                Ditot += Di
                Ytot += Y
                Ltot += L
        frc = self.nfres.nffrc.sum()
        mom = self.nfres.nfmom.sum()
        table1.add_row(['Total', frc.x, frc.y, frc.z, mom.x, mom.y, mom.z])
        table = report.add_table()
        table.add_column('Density', '.3f', data=[self.rho])
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Dynamic Pressure', '.1f', data=[self.qfs])
        table2.add_row(['Total', self.sys.sref, Ditot, Ytot, Ltot, self.nfres.CDi, self.nfres.CY, self.nfres.CL])
        return report

    @property
    def strip_forces(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('#', 'd')
        table.add_column('ypos', '.4f')
        table.add_column('Chord', '.4f')
        table.add_column('Area', '.4f')
        table.add_column('c cl', '.4f')
        table.add_column('ai', '.4f')
        table.add_column('cl_norm', '.4f')
        table.add_column('cl', '.4f')
        table.add_column('cd', '.4f')
        q = self.qfs
        for strp in self.sys.strps:
            j = strp.lsid
            ypos = strp.pnti.y
            chord = strp.chord
            area = strp.area
            nrmfrc = Vector(0.0, 0.0, 0.0)
            for pnl in strp.pnls:
                nrmfrc += self.nfres.nffrc[pnl.lpid]/pnl.area*pnl.crd
            c_cl = nrmfrc.z/q
            ai = -self.trres.trwsh[strp.lsid]/self.speed
            cd = self.trres.trfrc.x[strp.lsid]/q/area
            cy = nrmfrc.dot(self.acs.diry)/q/chord
            cz = nrmfrc.dot(self.acs.dirz)/q/chord
            cf = Vector(0.0, cy, cz)
            cl_norm = cf.dot(strp.nrmt)
            cl = cl_norm
            table.add_row([j, ypos, chord, area, c_cl, ai, cl_norm, cl, cd])
        return report

    @property
    def strip_coefficients(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('#', 'd')
        table.add_column('Chord', '.4f')
        table.add_column('Area', '.6f')
        table.add_column('cn', '.5f')
        table.add_column('ca', '.5f')
        table.add_column('cl', '.5f')
        table.add_column('cd', '.5f')
        table.add_column('dw', '.5f')
        table.add_column('cm le', '.5f')
        table.add_column('cm qc', '.5f')
        q = self.qfs
        for strp in self.sys.strps:
            j = strp.lsid
            chord = strp.chord
            area = strp.area
            force = Vector(0.0, 0.0, 0.0)
            momle = Vector(0.0, 0.0, 0.0)
            for pnl in strp.pnls:
                force += self.nfres.nffrc[pnl.lpid]
                rref = pnl.pnti - strp.pnti
                momle += rref.cross(self.nfres.nffrc[pnl.lpid])
            cn = force.dot(strp.nrmt)/q/area
            ca = force.x/q/area
            cl = force.dot(self.acs.dirz)/q/area
            cd = force.dot(self.acs.dirx)/q/area
            dw = -self.trres.trwsh[strp.lsid]/self.speed
            cmle = momle.y/q/area/chord
            rqc = Vector(-chord/4, 0.0, 0.0)
            momqc = momle + rqc.cross(force)
            cmqc = momqc.y/q/area/chord
            table.add_row([j, chord, area, cn, ca, cl, cd, dw, cmle, cmqc])
        return report

    @property
    def panel_forces(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('# P', 'd')
        table.add_column('# S', 'd')
        table.add_column('X', '.5f')
        table.add_column('Y', '.5f')
        table.add_column('Z', '.5f')
        table.add_column('DX', '.5f')
        table.add_column('Slope', '.5f')
        table.add_column('dCp', '.5f')
        q = self.qfs
        for pnl in self.sys.pnls:
            j = pnl.lpid
            k = pnl.strp.lsid
            x = pnl.pntg.x
            y = pnl.pntg.y
            z = pnl.pntg.z
            area = pnl.area
            frc = self.nfres.nffrc[pnl.lpid]
            nfrc = frc.dot(pnl.nrml)
            cp = nfrc/area/q
            chord = pnl.crd
            alc = tan(radians(pnl.alpha))
            table.add_row([j, k, x, y, z, chord, alc, cp])
        return report

    @property
    def panel_near_field_results(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('# P', 'd')
        table.add_column('# S', 'd')
        table.add_column('Gamma', '.1f')
        table.add_column('Vx', '.1f')
        table.add_column('Vy', '.1f')
        table.add_column('Vz', '.1f')
        table.add_column('lx', '.4f')
        table.add_column('ly', '.4f')
        table.add_column('lz', '.4f')
        table.add_column('Fx', '.2f')
        table.add_column('Fy', '.2f')
        table.add_column('Fz', '.2f')
        for pnl in self.sys.pnls:
            j = pnl.lpid
            k = pnl.strp.lsid
            gamma = self.gamma[j]
            Vx = self.nfres.nfvel[j].x
            Vy = self.nfres.nfvel[j].y
            Vz = self.nfres.nfvel[j].z
            lx = pnl.leni.x
            ly = pnl.leni.y
            lz = pnl.leni.z
            Fx = self.nfres.nffrc[j].x
            Fy = self.nfres.nffrc[j].y
            Fz = self.nfres.nffrc[j].z
            table.add_row([j, k, gamma, Vx, Vy, Vz, lx, ly, lz, Fx, Fy, Fz])
        return report

    @property
    def stability_derivatives(self) -> MDReport:
        return self.stres.stability_derivatives

    @property
    def stability_derivatives_body(self) -> MDReport:
        return self.stres.stability_derivatives_body

    @property
    def control_derivatives(self) -> MDReport:
        from . import sfrm
        report = MDReport()
        report.add_heading('Control Derivatives', 1)
        for control in self.ctrls:
            letter = control[0]
            report.add_heading(f'{control.capitalize()} Derivatives', 2)
            table = report.add_table()
            table.add_column(f'CLd{letter:s}', sfrm)
            table.add_column(f'CYd{letter:s}', sfrm)
            table.add_column(f'Cld{letter:s}', sfrm)
            table.add_column(f'Cmd{letter:s}', sfrm)
            table.add_column(f'Cnd{letter:s}', sfrm)
            if self.ctrls[control] >= 0.0:
                ctresp = self.ctresp[control]
                table.add_row([ctresp.CL, ctresp.CY, ctresp.Cl, ctresp.Cm, ctresp.Cn])
            if self.ctrls[control] <= 0.0:
                ctresp = self.ctresp[control]
                table.add_row([ctresp.CL, ctresp.CY, ctresp.Cl, ctresp.Cm, ctresp.Cn])
        return report

    def to_mdobj(self) -> MDReport:
        from . import cfrm, dfrm, efrm

        report = MDReport()
        report.add_heading(f'Lattice Result {self.name:s} for {self.sys.name:s}', 1)

        table = report.add_table()
        table.add_column('Alpha (deg)', cfrm, data=[self.alpha])
        table.add_column('Beta (deg)', cfrm, data=[self.beta])
        table.add_column('Speed', cfrm, data=[self.speed])
        table.add_column('Rho', cfrm, data=[self.rho])
        table.add_column('Mach', efrm, data=[self.mach])

        table = report.add_table()
        table.add_column('pb/2V (rad)', cfrm, data=[self.pbo2V])
        table.add_column('qc/2V (rad)', cfrm, data=[self.qco2V])
        table.add_column('rb/2V (rad)', cfrm, data=[self.rbo2V])

        table = report.add_table()
        table.add_column('xcg', '.5f', data=[self.rcg.x])
        table.add_column('ycg', '.5f', data=[self.rcg.y])
        table.add_column('zcg', '.5f', data=[self.rcg.z])

        if len(self.ctrls) > 0:
            table = report.add_table()
            for control in self.ctrls:
                ctrl = self.ctrls[control]
                control = control.capitalize()
                table.add_column(f'{control} (deg)', cfrm, data=[ctrl])

        if self.sys.cdo != 0.0:
            table = report.add_table()
            table.add_column('CDo', dfrm, data=[self.pdres.CDo])
            table.add_column('CYo', cfrm, data=[self.pdres.CY])
            table.add_column('CLo', cfrm, data=[self.pdres.CL])
            table.add_column('Clo', cfrm, data=[self.pdres.Cl])
            table.add_column('Cmo', cfrm, data=[self.pdres.Cm])
            table.add_column('Cno', cfrm, data=[self.pdres.Cn])

        if self.gamma is not None:
            table = report.add_table()
            table.add_column('Cx', cfrm, data=[self.nfres.Cx])
            table.add_column('Cy', cfrm, data=[self.nfres.Cy])
            table.add_column('Cz', cfrm, data=[self.nfres.Cz])

            table = report.add_table()
            table.add_column('CDi', dfrm, data=[self.nfres.CDi])
            table.add_column('CY', cfrm, data=[self.nfres.CY])
            table.add_column('CL', cfrm, data=[self.nfres.CL])
            table.add_column('Cl', cfrm, data=[self.nfres.Cl])
            table.add_column('Cm', cfrm, data=[self.nfres.Cm])
            table.add_column('Cn', cfrm, data=[self.nfres.Cn])

            if self.sys.cdo != 0.0:
                lod = self.nfres.CL/(self.pdres.CDo+self.nfres.CDi)
                table.add_column('L/D', '.5g', data=[lod])

        if self.phi is not None:
            table = report.add_table()
            table.add_column('CDi_ff', dfrm, data=[self.trres.CDi])
            table.add_column('CY_ff', cfrm, data=[self.trres.CY])
            table.add_column('CL_ff', cfrm, data=[self.trres.CL])
            table.add_column('e', efrm, data=[self.trres.e])
            if self.sys.cdo != 0.0:
                lod_ff = self.trres.CL/(self.pdres.CDo+self.trres.CDi)
                table.add_column('L/D_ff', '.5g', data=[lod_ff])

        return report

    @classmethod
    def from_dict(cls, sys: 'System', resdict: dict[str, Any]) -> 'LatticeResult':
        name = resdict.get('name', 'Lattice Result')
        inherit = resdict.get('inherit', None)
        if inherit is not None:
            result = sys.results[inherit].to_result(name=name)
        else:
            result = cls(name, sys)
        # Density
        rho = resdict.get('density', 1.0)
        result.set_density(rho=rho)
        # State
        mach = resdict.get('mach', None)
        speed = resdict.get('speed', None)
        alpha = resdict.get('alpha', None)
        beta = resdict.get('beta', None)
        pbo2V = resdict.get('pbo2V', None)
        qco2V = resdict.get('qco2V', None)
        rbo2V = resdict.get('rbo2V', None)
        result.set_state(mach=mach, speed=speed, alpha=alpha, beta=beta,
                         pbo2V=pbo2V, qco2V=qco2V, rbo2V=rbo2V)
        # Controls
        ctrls: dict[str, float] = {}
        for control in sys.ctrls:
            value = resdict.get(control, None)
            if value is not None:
                ctrls[control] = value
        result.set_controls(**ctrls)
        # Centre of Gravity
        rcgdict = resdict.get('rcg', None)
        if rcgdict is not None:
            rcg = Vector(rcgdict['x'], rcgdict['y'], rcgdict['z'])
            result.set_cg(rcg)
        sys.results[name] = result

        mass = resdict.get('mass', None)
        if isinstance(mass, dict):
            mass = Mass(**mass)
        elif isinstance(mass, float):
            mass = Mass(name = result.name, mass = mass)
        elif mass is None:
            if sys.mass is not None:
                mass = sys.mass
            else:
                mass = Mass(result.name, mass = 1.0, xcm = result.rcg.x,
                            ycm = result.rcg.y, zcm = result.rcg.z)
        result.mass = mass

        return result

    def __str__(self) -> str:
        return self.to_mdobj()._repr_markdown_()

    def __repr__(self) -> str:
        return f'<LatticeResult: {self.name}>'

    def _repr_markdown_(self) -> str:
        return self.to_mdobj()._repr_markdown_()


class GammaResult():
    res: LatticeResult = None
    gamma: 'NDArray' = None
    _rhogamma: 'NDArray' = None
    _nfvel: Vector = None
    _nffrc: Vector = None
    _nfmom: Vector = None
    _nffrctot: Vector = None
    _nfmomtot: Vector = None
    _Cx: float = None
    _Cy: float = None
    _Cz: float = None
    _Cmx: float = None
    _Cmy: float = None
    _Cmz: float = None
    _CDi: float = None
    _CY: float = None
    _CL: float = None
    _e: float = None
    _Cl: float = None
    _Cm: float = None
    _Cn: float = None

    def __init__(self, res: LatticeResult, gamma: 'NDArray') -> None:
        self.res = res
        self.gamma = gamma

    @property
    def rhogamma(self) -> 'NDArray':
        if self._rhogamma is None:
            self._rhogamma = self.res.rho*self.gamma
        return self._rhogamma

    @property
    def nfvel(self) -> Vector:
        if self._nfvel is None:
            self._nfvel = self.res.avg@self.gamma + self.res.avv
        return self._nfvel

    @property
    def nffrc(self) -> Vector:
        if self._nffrc is None:
            tmp = self.res.afg@self.gamma + self.res.afv
            self._nffrc = tmp*self.rhogamma
        return self._nffrc

    @property
    def nfmom(self) -> Vector:
        if self._nfmom is None:
            self._nfmom = self.res.arm.cross(self.nffrc)
        return self._nfmom

    @property
    def nffrctot(self) -> Vector:
        if self._nffrctot is None:
            self._nffrctot = self.nffrc.sum()
        return self._nffrctot

    @property
    def nfmomtot(self) -> Vector:
        if self._nfmomtot is None:
            self._nfmomtot = self.nfmom.sum()
        return self._nfmomtot

    @property
    def Cx(self) -> float:
        if self._Cx is None:
            self._Cx = self.nffrctot.x/self.res.qfs/self.res.sys.sref
            self._Cx = fix_zero(self._Cx)
        return self._Cx

    @property
    def Cy(self) -> float:
        if self._Cy is None:
            self._Cy = self.nffrctot.y/self.res.qfs/self.res.sys.sref
            self._Cy = fix_zero(self._Cy)
        return self._Cy

    @property
    def Cz(self) -> float:
        if self._Cz is None:
            self._Cz = self.nffrctot.z/self.res.qfs/self.res.sys.sref
            self._Cz = fix_zero(self._Cz)
        return self._Cz

    @property
    def Cmx(self) -> float:
        if self._Cmx is None:
            self._Cmx = self.nfmomtot.x/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cmx = fix_zero(self._Cmx)
        return self._Cmx

    @property
    def Cmy(self) -> float:
        if self._Cmy is None:
            self._Cmy = self.nfmomtot.y/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cmy = fix_zero(self._Cmy)
        return self._Cmy

    @property
    def Cmz(self) -> float:
        if self._Cmz is None:
            self._Cmz = self.nfmomtot.z/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cmz = fix_zero(self._Cmz)
        return self._Cmz

    @property
    def CDi(self) -> float:
        if self._CDi is None:
            Di = self.nffrctot.dot(self.res.acs.dirx)
            self._CDi = Di/self.res.qfs/self.res.sys.sref
            self._CDi = fix_zero(self._CDi)
        return self._CDi

    @property
    def CY(self) -> float:
        if self._CY is None:
            Y = self.nffrctot.dot(self.res.acs.diry)
            self._CY = Y/self.res.qfs/self.res.sys.sref
            self._CY = fix_zero(self._CY)
        return self._CY

    @property
    def CL(self) -> float:
        if self._CL is None:
            L = self.nffrctot.dot(self.res.acs.dirz)
            self._CL = L/self.res.qfs/self.res.sys.sref
            self._CL = fix_zero(self._CL)
        return self._CL

    @property
    def Cl(self) -> float:
        if self._Cl is None:
            l = self.nfmomtot.dot(self.res.scs.dirx)
            self._Cl = l/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cl = fix_zero(self._Cl)
        return self._Cl

    @property
    def Cm(self) -> float:
        if self._Cm is None:
            m = self.nfmomtot.dot(self.res.scs.diry)
            self._Cm = m/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cm = fix_zero(self._Cm)
        return self._Cm

    @property
    def Cn(self) -> float:
        if self._Cn is None:
            n = self.nfmomtot.dot(self.res.scs.dirz)
            self._Cn = n/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cn = fix_zero(self._Cn)
        return self._Cn

    @property
    def e(self) -> float:
        if self._e is None:
            if self.CDi <= 0.0:
                self._e = float('nan')
            elif self.CL == 0.0 and self.CY == 0.0:
                self._e = 0.0
            else:
                self._e = (self.CL**2+self.CY**2)/pi/self.res.sys.ar/self.CDi
                self._e = fix_zero(self._e)
        return self._e


class StabilityGammaResult:
    res: LatticeResult = None
    dvfs: Vector = None
    dofs: Vector = None
    dacs: dict[str, Vector] = None
    dscs: dict[str, Vector] = None
    _dgamma: 'NDArray' = None
    _rhodgamma: 'NDArray' = None
    _davv: Vector = None
    _dafv: Vector = None
    _dvel: Vector = None
    _dfrc: Vector = None
    _dmom: Vector = None
    _dfrctot: Vector = None
    _dmomtot: Vector = None
    _Cx: float = None
    _Cy: float = None
    _Cz: float = None
    _Cmx: float = None
    _Cmy: float = None
    _Cmz: float = None
    _CDi: float = None
    _CY: float = None
    _CL: float = None
    _e: float = None
    _Cl: float = None
    _Cm: float = None
    _Cn: float = None

    def __init__(self, res: LatticeResult, dvfs: Vector | None = None, dofs: Vector | None = None,
                 dacs: dict[str, Vector] = None, dscs: dict[str, Vector] | None = None) -> None:
        self.res = res
        self.dvfs = dvfs
        self.dofs = dofs
        self.dacs = dacs
        self.dscs = dscs

    @property
    def dgamma(self) -> 'NDArray':
        if self._dgamma is None:
            self._dgamma = zeros(self.res.gamma.shape)
            if self.dvfs is not None:
                self._dgamma += self.res.ungam[:, 0].dot(self.dvfs)
            if self.dofs is not None:
                self._dgamma += self.res.ungam[:, 1].dot(self.dofs)
        return self._dgamma

    @property
    def rhodgamma(self) -> 'NDArray':
        if self._rhodgamma is None:
            self._rhodgamma = self.res.rho*self.dgamma
        return self._rhodgamma

    @property
    def davv(self) -> Vector:
        if self._davv is None:
            num = len(self.res.sys.pnls)
            self._davv = Vector.zeros(num)
            for pnl in self.res.sys.pnls:
                i = pnl.lpid
                if self.dvfs is not None:
                    self._davv[i] += self.dvfs
                if self.dofs is not None:
                    self._davv[i] -= self.dofs.cross(self.res.arm[i])
        return self._davv

    @property
    def dafv(self) -> Vector:
        if self._dafv is None:
            num = len(self.res.sys.pnls)
            self._dafv = Vector.zeros(num)
            for pnl in self.res.sys.pnls:
                if not pnl.noload:
                    i = pnl.lpid
                    self._dafv[i] = self.davv[i].cross(pnl.leni)
        return self._dafv

    @property
    def dvel(self) -> Vector:
        if self._dvel is None:
            self._dvel = self.res.avg@self.dgamma + self.davv
        return self._dvel

    @property
    def dfrc(self) -> Vector:
        if self._dfrc is None:
            self._dfrc = (self.res.afg@self.dgamma + self.dafv)*self.res.nfres.rhogamma
            self._dfrc += (self.res.afg@self.res.gamma + self.res.afv)*self.rhodgamma
        return self._dfrc

    @property
    def dmom(self) -> Vector:
        if self._dmom is None:
            self._dmom = self.res.arm.cross(self.dfrc)
        return self._dmom

    @property
    def dfrctot(self) -> Vector:
        if self._dfrctot is None:
            self._dfrctot = self.dfrc.sum()
        return self._dfrctot

    @property
    def dmomtot(self) -> Vector:
        if self._dmomtot is None:
            self._dmomtot = self.dmom.sum()
        return self._dmomtot

    @property
    def Cx(self) -> float:
        if self._Cx is None:
            self._Cx = self.dfrctot.x/self.res.qfs/self.res.sys.sref
            self._Cx = fix_zero(self._Cx)
        return self._Cx

    @property
    def Cy(self) -> float:
        if self._Cy is None:
            self._Cy = self.dfrctot.y/self.res.qfs/self.res.sys.sref
            self._Cy = fix_zero(self._Cy)
        return self._Cy

    @property
    def Cz(self) -> float:
        if self._Cz is None:
            self._Cz = self.dfrctot.z/self.res.qfs/self.res.sys.sref
            self._Cz = fix_zero(self._Cz)
        return self._Cz

    @property
    def Cmx(self) -> float:
        if self._Cmx is None:
            self._Cmx = self.dmomtot.x/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cmx = fix_zero(self._Cmx)
        return self._Cmx

    @property
    def Cmy(self) -> float:
        if self._Cmy is None:
            self._Cmy = self.dmomtot.y/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cmy = fix_zero(self._Cmy)
        return self._Cmy

    @property
    def Cmz(self) -> float:
        if self._Cmz is None:
            self._Cmz = self.dmomtot.z/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cmz = fix_zero(self._Cmz)
        return self._Cmz

    @property
    def CDi(self) -> float:
        if self._CDi is None:
            Di = self.dfrctot.dot(self.res.acs.dirx)
            if self.dacs is not None:
                Di += self.res.nfres.nffrctot.dot(self.dacs['x'])
            self._CDi = Di/self.res.qfs/self.res.sys.sref
            self._CDi = fix_zero(self._CDi)
        return self._CDi

    @property
    def CY(self) -> float:
        if self._CY is None:
            Y = self.dfrctot.dot(self.res.acs.diry)
            if self.dacs is not None:
                Y += self.res.nfres.nffrctot.dot(self.dacs['y'])
            self._CY = Y/self.res.qfs/self.res.sys.sref
            self._CY = fix_zero(self._CY)
        return self._CY

    @property
    def CL(self) -> float:
        if self._CL is None:
            L = self.dfrctot.dot(self.res.acs.dirz)
            if self.dacs is not None:
                L += self.res.nfres.nffrctot.dot(self.dacs['z'])
            self._CL = L/self.res.qfs/self.res.sys.sref
            self._CL = fix_zero(self._CL)
        return self._CL

    @property
    def Cl(self) -> float:
        if self._Cl is None:
            l = self.dmomtot.dot(self.res.scs.dirx)
            if self.dscs is not None:
                l += self.res.nfres.nfmomtot.dot(self.dscs['x'])
            self._Cl = l/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cl = fix_zero(self._Cl)
        return self._Cl

    @property
    def Cm(self) -> float:
        if self._Cm is None:
            m = self.dmomtot.dot(self.res.scs.diry)
            if self.dscs is not None:
                m += self.res.nfres.nfmomtot.dot(self.dscs['y'])
            self._Cm = m/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cm = fix_zero(self._Cm)
        return self._Cm

    @property
    def Cn(self) -> float:
        if self._Cn is None:
            n = self.dmomtot.dot(self.res.scs.dirz)
            if self.dscs is not None:
                n += self.res.nfres.nfmomtot.dot(self.dscs['z'])
            self._Cn = n/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cn = fix_zero(self._Cn)
        return self._Cn


class StripResult():
    gamres: GammaResult = None
    _stfrc: Vector = None
    _stmom: Vector = None
    _lift: float = None
    _side: float = None
    _drag: float = None
    _pmom: float = None

    def __init__(self, gamres: GammaResult) -> None:
        self.gamres = gamres

    @property
    def stfrc(self) -> Vector:
        if self._stfrc is None:
            sys = self.gamres.res.sys
            num = len(sys.strps)
            self._stfrc = Vector.zeros(num)
            for strp in sys.strps:
                i = strp.lsid
                for pnl in strp.pnls:
                    j = pnl.lpid
                    self._stfrc[i] += self.gamres.nffrc[j]
        return self._stfrc

    @property
    def stmom(self) -> Vector:
        if self._stmom is None:
            sys = self.gamres.res.sys
            num = len(sys.strps)
            self._stmom = Vector.zeros(num)
            for strp in sys.strps:
                i = strp.lsid
                for pnl in strp.pnls:
                    j = pnl.lpid
                    arm = pnl.pnti - strp.pntq
                    self._stmom[i] += arm.cross(self.gamres.nffrc[j])
        return self._stmom

    @property
    def drag(self) -> float:
        if self._drag is None:
            self._drag = self.stfrc.dot(self.gamres.res.acs.dirx)
        return self._drag

    @property
    def side(self) -> float:
        if self._side is None:
            self._side = self.stfrc.dot(self.gamres.res.acs.diry)
        return self._side

    @property
    def lift(self) -> float:
        if self._lift is None:
            self._lift = self.stfrc.dot(self.gamres.res.acs.dirz)
        return self._lift

    @property
    def pmom(self) -> float:
        if self._pmom is None:
            self._pmom = self.stmom.dot(self.gamres.res.acs.diry)
        return self._pmom


class PhiResult():
    res: LatticeResult = None
    phi: 'NDArray' = None
    _trwsh: 'NDArray' = None
    _trfrc: Vector = None
    _trmom: Vector = None
    _trfrctot: Vector = None
    _trmomtot: Vector = None
    _CDi: float = None
    _CY: float = None
    _CL: float = None
    _Cl: float = None
    _Cm: float = None
    _Cn: float = None
    _e: float = None
    _lod: float = None

    def __init__(self, res: LatticeResult, phi: 'NDArray') -> None:
        self.res = res
        self.phi = phi

    @property
    def trwsh(self):
        if self._trwsh is None:
            self._trwsh = self.res.sys.bvg@self.phi
        return self._trwsh

    @property
    def trfrc(self):
        if self._trfrc is None:
            x = self.res.rho*self.phi*(self.res.sys.bdg@self.phi)
            y = self.res.rho*self.res.speed*self.phi*self.res.sys.byg
            z = self.res.rho*self.res.speed*self.phi*self.res.sys.blg
            self._trfrc = Vector(x, y, z)
        return self._trfrc

    @property
    def trmom(self):
        if self._trmom is None:
            self._trmom = self.res.brm.cross(self.trfrc)
        return self._trmom

    @property
    def trfrctot(self) -> Vector:
        if self._trfrctot is None:
            self._trfrctot = self.trfrc.sum()
        return self._trfrctot

    @property
    def trmomtot(self) -> Vector:
        if self._trmomtot is None:
            self._trmomtot = self.trmom.sum()
        return self._trmomtot

    @property
    def CDi(self) -> float:
        if self._CDi is None:
            Di = self.trfrctot.x
            self._CDi = Di/self.res.qfs/self.res.sys.sref
            self._CDi = fix_zero(self._CDi)
        return self._CDi

    @property
    def CY(self) -> float:
        if self._CY is None:
            Y = self.trfrctot.y
            self._CY = Y/self.res.qfs/self.res.sys.sref
            self._CY = fix_zero(self._CY)
        return self._CY

    @property
    def CL(self) -> float:
        if self._CL is None:
            L = self.trfrctot.z
            self._CL = L/self.res.qfs/self.res.sys.sref
            self._CL = fix_zero(self._CL)
        return self._CL

    @property
    def Cl(self) -> float:
        if self._Cl is None:
            l = -self.trmomtot.x
            self._Cl = l/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cl = fix_zero(self._Cl)
        return self._Cl

    @property
    def Cm(self) -> float:
        if self._Cm is None:
            m = self.trmomtot.y
            self._Cm = m/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cm = fix_zero(self._Cm)
        return self._Cm

    @property
    def Cn(self) -> float:
        if self._Cn is None:
            n = -self.trmomtot.z
            self._Cn = n/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cn = fix_zero(self._Cn)
        return self._Cn

    @property
    def e(self) -> float:
        if self._e is None:
            if self.CDi == 0.0:
                if self.CL == 0.0 and self.CY == 0.0:
                    self._e = 0.0
                else:
                    self._e = float('nan')
            else:
                self._e = (self.CL**2+self.CY**2)/pi/self.res.sys.ar/self.CDi
                self._e = fix_zero(self._e)
        return self._e


class ParasiticDragResult():
    res: LatticeResult = None
    _pdfrc: Vector = None
    _pdmom: Vector = None
    _pdfrctot: Vector = None
    _pdmomtot: Vector = None
    _CDo: float = None
    _CY: float = None
    _CL: float = None
    _Cl: float = None
    _Cm: float = None
    _Cn: float = None

    def __init__(self, res: LatticeResult):
        self.res = res

    @property
    def pdfrc(self) -> Vector:
        if self._pdfrc is None:
            dynpr = self.res.bvv.dot(self.res.bvv)*(self.res.rho/2)
            tmp = dynpr*self.res.sys.bda
            self._pdfrc = Vector(tmp*self.res.acs.dirx.x,
                                      tmp*self.res.acs.dirx.y,
                                      tmp*self.res.acs.dirx.z)
        return self._pdfrc

    @property
    def pdmom(self) -> Vector:
        if self._pdmom is None:
            self._pdmom = self.res.brm.cross(self.pdfrc)
        return self._pdmom

    @property
    def pdfrctot(self) -> Vector:
        if self._pdfrctot is None:
            self._pdfrctot = self.pdfrc.sum()
        return self._pdfrctot

    @property
    def pdmomtot(self) -> Vector:
        if self._pdmomtot is None:
            self._pdmomtot = self.pdmom.sum()
        return self._pdmomtot

    @property
    def CDo(self) -> float:
        if self._CDo is None:
            Do = self.res.acs.dirx.dot(self.pdfrctot)
            self._CDo = Do/self.res.qfs/self.res.sys.sref
            self._CDo = fix_zero(self._CDo)
        return self._CDo

    @property
    def CY(self) -> float:
        if self._CY is None:
            Y = self.res.acs.diry.dot(self.pdfrctot)
            self._CY = Y/self.res.qfs/self.res.sys.sref
            self._CY = fix_zero(self._CY)
        return self._CY

    @property
    def CL(self) -> float:
        if self._CL is None:
            L = self.res.acs.dirz.dot(self.pdfrctot)
            self._CL = L/self.res.qfs/self.res.sys.sref
            self._CL = fix_zero(self._CL)
        return self._CL

    @property
    def Cl(self) -> float:
        if self._Cl is None:
            l = self.res.scs.dirx.dot(self.pdmomtot)
            self._Cl = l/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cl = fix_zero(self._Cl)
        return self._Cl

    @property
    def Cm(self) -> float:
        if self._Cm is None:
            m = self.res.scs.diry.dot(self.pdmomtot)
            self._Cm = m/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cm = fix_zero(self._Cm)
        return self._Cm

    @property
    def Cn(self) -> float:
        if self._Cn is None:
            n = self.res.scs.dirz.dot(self.pdmomtot)
            self._Cn = n/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cn = fix_zero(self._Cn)
        return self._Cn


class StabilityResult():
    res: LatticeResult = None
    _u: StabilityGammaResult = None
    _v: StabilityGammaResult = None
    _w: StabilityGammaResult = None
    _p: StabilityGammaResult = None
    _q: StabilityGammaResult = None
    _r: StabilityGammaResult = None
    _alpha: StabilityGammaResult = None
    _beta: StabilityGammaResult = None
    _pbo2V: StabilityGammaResult = None
    _qco2V: StabilityGammaResult = None
    _rbo2V: StabilityGammaResult = None
    _pdbo2V: StabilityGammaResult = None
    _qdco2V: StabilityGammaResult = None
    _rdbo2V: StabilityGammaResult = None
    _xnp: float = None
    _sprat: float = None

    def __init__(self, res: LatticeResult) -> None:
        self.res = res

    @property
    def u(self) -> StabilityGammaResult:
        if self._u is None:
            dvfs = Vector(1.0, 0.0, 0.0)
            self._u = StabilityGammaResult(self.res, dvfs = dvfs)
        return self._u

    @property
    def v(self) -> StabilityGammaResult:
        if self._v is None:
            dvfs = Vector(0.0, 1.0, 0.0)
            self._v = StabilityGammaResult(self.res, dvfs = dvfs)
        return self._v

    @property
    def w(self) -> StabilityGammaResult:
        if self._w is None:
            dvfs = Vector(0.0, 0.0, 1.0)
            self._w = StabilityGammaResult(self.res, dvfs = dvfs)
        return self._w

    @property
    def p(self) -> StabilityGammaResult:
        if self._p is None:
            dpqr = Vector(1.0, 0.0, 0.0)
            dofs = self.res.scs.vector_to_global(dpqr)
            self._p = StabilityGammaResult(self.res, dofs = dofs)
        return self._p

    @property
    def q(self) -> StabilityGammaResult:
        if self._q is None:
            dpqr = Vector(0.0, 1.0, 0.0)
            dofs = self.res.scs.vector_to_global(dpqr)
            self._q = StabilityGammaResult(self.res, dofs = dofs)
        return self._q

    @property
    def r(self) -> StabilityGammaResult:
        if self._r is None:
            dpqr = Vector(0.0, 0.0, 1.0)
            dofs = self.res.scs.vector_to_global(dpqr)
            self._r = StabilityGammaResult(self.res, dofs = dofs)
        return self._r

    @property
    def alpha(self) -> StabilityGammaResult:
        if self._alpha is None:
            dvfs = self.res.dvfsa
            dofs = Vector(
                self.res.dscsa['x'].dot(self.res.pqr),
                self.res.dscsa['y'].dot(self.res.pqr),
                self.res.dscsa['z'].dot(self.res.pqr)
            )
            self._alpha = StabilityGammaResult(self.res, dvfs = dvfs, dofs = dofs,
                                               dacs = self.res.dacsa,
                                               dscs = self.res.dscsa)
        return self._alpha

    @property
    def beta(self) -> StabilityGammaResult:
        if self._beta is None:
            dvfs = self.res.dvfsb
            self._beta = StabilityGammaResult(self.res, dvfs = dvfs,
                                              dacs = self.res.dacsb,
                                              dscs = self.res.dscsb)
        return self._beta

    @property
    def pbo2V(self) -> StabilityGammaResult:
        if self._pbo2V is None:
            dpqr = Vector(2*self.res.speed/self.res.sys.bref, 0.0, 0.0)
            dofs = self.res.scs.vector_to_global(dpqr)
            self._pbo2V = StabilityGammaResult(self.res, dofs = dofs)
        return self._pbo2V

    @property
    def qco2V(self) -> StabilityGammaResult:
        if self._qco2V is None:
            dpqr = Vector(0.0, 2*self.res.speed/self.res.sys.cref, 0.0)
            dofs = self.res.scs.vector_to_global(dpqr)
            self._qco2V = StabilityGammaResult(self.res, dofs = dofs)
        return self._qco2V

    @property
    def rbo2V(self) -> StabilityGammaResult:
        if self._rbo2V is None:
            dpqr = Vector(0.0, 0.0, 2*self.res.speed/self.res.sys.bref)
            dofs = self.res.scs.vector_to_global(dpqr)
            self._rbo2V = StabilityGammaResult(self.res, dofs = dofs)
        return self._rbo2V

    @property
    def pdbo2V(self) -> StabilityGammaResult:
        if self._pdbo2V is None:
            dofs = Vector(2*self.res.speed/self.res.sys.bref, 0.0, 0.0)
            self._pdbo2V = StabilityGammaResult(self.res, dofs = dofs)
        return self._pdbo2V

    @property
    def qdco2V(self) -> GammaResult:
        if self._qdco2V is None:
            dofs = Vector(0.0, 2*self.res.speed/self.res.sys.cref, 0.0)
            self._qdco2V = StabilityGammaResult(self.res, dofs = dofs)
        return self._qdco2V

    @property
    def rdbo2V(self) -> GammaResult:
        if self._rdbo2V is None:
            dofs = Vector(0.0, 0.0, 2*self.res.speed/self.res.sys.bref)
            self._rdbo2V = StabilityGammaResult(self.res, dofs = dofs)
        return self._rdbo2V

    @property
    def xnp(self) -> float:
        if self._xnp is None:
            xcg = self.res.rcg.x
            CLa = self.alpha.CL
            CMa = self.alpha.Cm
            c = self.res.sys.cref
            self._xnp = xcg - c*CMa/CLa
        return self._xnp

    @property
    def sprat(self) -> float:
        if self._sprat is None:
            Clb = self.beta.Cl
            Cnb = self.beta.Cn
            Cnr = self.rbo2V.Cn
            Clr = self.rbo2V.Cl
            if Clb == 0.0 and Clr == 0.0:
                self._sprat = float('nan')
            elif Cnb == 0.0 and Cnr == 0.0:
                self._sprat = float('nan')
            else:
                self._sprat = Clb*Cnr/(Clr*Cnb)
        return self._sprat

    def system_aerodynamic_matrix(self) -> 'NDArray':
        A = zeros((6, 6))
        F = self.u.dfrctot
        A[0, 0], A[1, 0], A[2, 0] = F.x, F.y, F.z
        M = self.res.scs.vector_to_local(self.u.dmomtot)
        A[3, 0], A[4, 0], A[5, 0] = M.x, M.y, M.z
        F = self.v.dfrctot
        A[0, 1], A[1, 1], A[2, 1] = F.x, F.y, F.z
        M = self.res.scs.vector_to_local(self.v.dmomtot)
        A[3, 1], A[4, 1], A[5, 1] = M.x, M.y, M.z
        F = self.w.dfrctot
        A[0, 2], A[1, 2], A[2, 2] = F.x, F.y, F.z
        M = self.res.scs.vector_to_local(self.w.dmomtot)
        A[3, 2], A[4, 2], A[5, 2] = M.x, M.y, M.z
        F = self.p.dfrctot
        A[0, 3], A[1, 3], A[2, 3] = F.x, F.y, F.z
        M = self.res.scs.vector_to_local(self.p.dmomtot)
        A[3, 3], A[4, 3], A[5, 3] = M.x, M.y, M.z
        F = self.q.dfrctot
        A[0, 4], A[1, 4], A[2, 4] = F.x, F.y, F.z
        M = self.res.scs.vector_to_local(self.q.dmomtot)
        A[3, 4], A[4, 4], A[5, 4] = M.x, M.y, M.z
        F = self.r.dfrctot
        A[0, 5], A[1, 5], A[2, 5] = F.x, F.y, F.z
        M = self.res.scs.vector_to_local(self.r.dmomtot)
        A[3, 5], A[4, 5], A[5, 5] = M.x, M.y, M.z
        return A

    @property
    def stability_derivatives(self) -> MDReport:

        from . import sfrm
        report = MDReport()
        report.add_heading('Stability Derivatives', 2)
        table = report.add_table()
        table.add_column('CLa', sfrm, data=[self.alpha.CL])
        table.add_column('CYa', sfrm, data=[self.alpha.CY])
        table.add_column('Cla', sfrm, data=[self.alpha.Cl])
        table.add_column('Cma', sfrm, data=[self.alpha.Cm])
        table.add_column('Cna', sfrm, data=[self.alpha.Cn])
        table = report.add_table()
        table.add_column('CLb', sfrm, data=[self.beta.CL])
        table.add_column('CYb', sfrm, data=[self.beta.CY])
        table.add_column('Clb', sfrm, data=[self.beta.Cl])
        table.add_column('Cmb', sfrm, data=[self.beta.Cm])
        table.add_column('Cnb', sfrm, data=[self.beta.Cn])
        table = report.add_table()
        table.add_column('CLp', sfrm, data=[self.pbo2V.CL])
        table.add_column('CYp', sfrm, data=[self.pbo2V.CY])
        table.add_column('Clp', sfrm, data=[self.pbo2V.Cl])
        table.add_column('Cmp', sfrm, data=[self.pbo2V.Cm])
        table.add_column('Cnp', sfrm, data=[self.pbo2V.Cn])
        table = report.add_table()
        table.add_column('CLq', sfrm, data=[self.qco2V.CL])
        table.add_column('CYq', sfrm, data=[self.qco2V.CY])
        table.add_column('Clq', sfrm, data=[self.qco2V.Cl])
        table.add_column('Cmq', sfrm, data=[self.qco2V.Cm])
        table.add_column('Cnq', sfrm, data=[self.qco2V.Cn])
        table = report.add_table()
        table.add_column('CLr', sfrm, data=[self.rbo2V.CL])
        table.add_column('CYr', sfrm, data=[self.rbo2V.CY])
        table.add_column('Clr', sfrm, data=[self.rbo2V.Cl])
        table.add_column('Cmr', sfrm, data=[self.rbo2V.Cm])
        table.add_column('Cnr', sfrm, data=[self.rbo2V.Cn])
        report.add_heading(f'Neutral Point Xnp = {self.xnp:.6f}', 3)
        report.add_heading(f'Clb.Cnr/(Clr.Cnb) = {self.sprat:.6f} (> 1 if spirally stable)', 3)
        return report

    @property
    def stability_derivatives_body(self) -> MDReport:
        from . import sfrm
        report = MDReport()
        report.add_heading('Stability Derivatives Body Axis', 2)
        table = report.add_table()
        table.add_column('Cxu', sfrm, data=[self.u.Cx])
        table.add_column('Cyu', sfrm, data=[self.u.Cy])
        table.add_column('Czu', sfrm, data=[self.u.Cz])
        table.add_column('Clu', sfrm, data=[self.u.Cmx])
        table.add_column('Cmu', sfrm, data=[self.u.Cmy])
        table.add_column('Cnu', sfrm, data=[self.u.Cmz])
        table = report.add_table()
        table.add_column('Cxv', sfrm, data=[self.v.Cx])
        table.add_column('Cyv', sfrm, data=[self.v.Cy])
        table.add_column('Czv', sfrm, data=[self.v.Cz])
        table.add_column('Clv', sfrm, data=[self.v.Cmx])
        table.add_column('Cmv', sfrm, data=[self.v.Cmy])
        table.add_column('Cnv', sfrm, data=[self.v.Cmz])
        table = report.add_table()
        table.add_column('Cxw', sfrm, data=[self.w.Cx])
        table.add_column('Cyw', sfrm, data=[self.w.Cy])
        table.add_column('Czw', sfrm, data=[self.w.Cz])
        table.add_column('Clw', sfrm, data=[self.w.Cmx])
        table.add_column('Cmw', sfrm, data=[self.w.Cmy])
        table.add_column('Cnw', sfrm, data=[self.w.Cmz])
        table = report.add_table()
        table.add_column('Cxp', sfrm, data=[self.pdbo2V.Cx])
        table.add_column('Cyp', sfrm, data=[self.pdbo2V.Cy])
        table.add_column('Czp', sfrm, data=[self.pdbo2V.Cz])
        table.add_column('Clp', sfrm, data=[self.pdbo2V.Cmx])
        table.add_column('Cmp', sfrm, data=[self.pdbo2V.Cmy])
        table.add_column('Cnp', sfrm, data=[self.pdbo2V.Cmz])
        table = report.add_table()
        table.add_column('Cxq', sfrm, data=[self.qdco2V.Cx])
        table.add_column('Cyq', sfrm, data=[self.qdco2V.Cy])
        table.add_column('Czq', sfrm, data=[self.qdco2V.Cz])
        table.add_column('Clq', sfrm, data=[self.qdco2V.Cmx])
        table.add_column('Cmq', sfrm, data=[self.qdco2V.Cmy])
        table.add_column('Cnq', sfrm, data=[self.qdco2V.Cmz])
        table = report.add_table()
        table.add_column('Cxr', sfrm, data=[self.rdbo2V.Cx])
        table.add_column('Cyr', sfrm, data=[self.rdbo2V.Cy])
        table.add_column('Czr', sfrm, data=[self.rdbo2V.Cz])
        table.add_column('Clr', sfrm, data=[self.rdbo2V.Cmx])
        table.add_column('Cmr', sfrm, data=[self.rdbo2V.Cmy])
        table.add_column('Cnr', sfrm, data=[self.rdbo2V.Cmz])
        report.add_object(table)
        return report

    def __str__(self) -> str:
        return self.stability_derivatives._repr_markdown_()

    def _repr_markdown_(self) -> str:
        return self.__str__()

def fix_zero(value: float, tol: float=1e-8):
    if abs(value) < tol:
        value = 0.0
    return value

def latticeresult_from_dict(lsys: 'System', resdata: dict[str, Any]) -> LatticeResult:
    result = LatticeResult.from_dict(lsys, resdata)
    return result

def trig_angle(angle: float) -> float:
    '''Calculates cos(angle) and sin(angle) with angle in degrees.'''
    angrad = radians(angle)
    cosang = cos(angrad)
    sinang = sin(angrad)
    return cosang, sinang
