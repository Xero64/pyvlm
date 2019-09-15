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
        # nrml = [strpi.leni.y for strpi in self.sys.strps]
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
            rho = self.rho
            tmp = self.sys.afg*self.gmat+self.afv
            self._nffrc = [gam*tmp[i, 0]*rho for i, gam in enumerate(self.gam)]
        return self._nffrc
    @property
    def nfmom(self):
        if self._nfmom is None:
            rho = self.rho
            tmp = self.sys.amg*self.gmat+self.amv
            self._nfmom = [gam*tmp[i, 0]*rho for i, gam in enumerate(self.gam)]
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
    def plot_panel_near_field_velocities(self, ax=None, component=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        py = [pnl.pnti.y for pnl in self.sys.pnls]
        if component is None or component == 'x':
            vx = [vel.x for vel in self.nfvel]
            ax.plot(py, vx, label=self.name+' Velocity X')
        if component is None or component == 'y':
            vy = [vel.y for vel in self.nfvel]
            ax.plot(py, vy, label=self.name+' Velocity Y')
        if component is None or component == 'z':
            vz = [vel.z for vel in self.nfvel]
            ax.plot(py, vz, label=self.name+' Velocity Z')
        ax.legend()
        return ax
    def plot_phi_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.sys.strpy, self.phi, label=self.name)
        ax.legend()
        return ax
    def plot_trefftz_lift_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        trlft = [self.trlft[i]/strp.dyt for i, strp in enumerate(self.sys.strps)]
        ax.plot(self.sys.strpy, trlft, label=self.name)
        ax.legend()
        return ax
    def plot_trefftz_drag_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        trdrg = [self.trdrg[i]/strp.dst for i, strp in enumerate(self.sys.strps)]
        ax.plot(self.sys.strpy, trdrg, label=self.name)
        ax.legend()
        return ax
    def plot_trefftz_wash_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.sys.strpy, self.trwsh, label=self.name)
        ax.legend()
        return ax
    def plot_strip_wash_distribution(self, lsid: int, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        stwash = []
        for i in range(self.sys.bvg.shape[0]):
            stwash.append(self.sys.bvg[i, lsid])
        ax.plot(self.sys.strpy, stwash)
        return ax
    def print_near_field_total_loads(self):
        print(f'Total Force in X = {self.nffrctot.x}')
        print(f'Total Force in Y = {self.nffrctot.y}')
        print(f'Total Force in Z = {self.nffrctot.z}')
        print(f'Total Moment in X = {self.nfmomtot.x}')
        print(f'Total Moment in Y = {self.nfmomtot.y}')
        print(f'Total Moment in Z = {self.nfmomtot.z}')
    def print_aerodynamic_coefficients(self):
        from . import cfrm, dfrm, efrm
        if self.gam is not None:
            print('Cx = '+self.Cx.__format__(cfrm))
            print('Cy = '+self.Cy.__format__(cfrm))
            print('Cz = '+self.Cz.__format__(cfrm))
            print('Cl = '+self.Cl.__format__(cfrm))
            print('Cm = '+self.Cm.__format__(cfrm))
            print('Cn = '+self.Cn.__format__(cfrm))
            print('CL = '+self.CL.__format__(cfrm))
            print('CDi = '+self.CDi.__format__(dfrm))
            print('CY = '+self.CY.__format__(cfrm))
        if self.phi is not None:
            print('CL_ff = '+self.CL_ff.__format__(cfrm))
            print('CDi_ff = '+self.CDi_ff.__format__(dfrm))
            print('CY_ff = '+self.CY_ff.__format__(cfrm))
            print('e = '+self.e.__format__(efrm))
    def print_strip_forces(self):
        from py2md.classes import MDTable
        from math import atan
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('Yle', '.4f')
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
            yle = strp.pnti.y
            chord = strp.chord
            area = strp.area
            nrmfrc = Vector(0.0, 0.0, 0.0)
            for pnl in strp.pnls:
                nrmfrc += self.nffrc[pnl.lpid]/pnl.area*pnl.crd
            c_cl = nrmfrc.z/q
            ai = -self.trwsh[strp.lsid]
            cd = self.trdrg[strp.lsid]/q/area
            # cx = nrmfrc*self.udc/q/chord
            cy = nrmfrc*self.uyc/q/chord
            cz = nrmfrc*self.ulc/q/chord
            cf = Vector(0.0, cy, cz)
            cl_norm = cf*strp.nrmt
            cl = cl_norm
            # cd = nrmfrc*self.udc/q/chord
            # cd = cx
            table.add_row([j, yle, chord, area, c_cl, ai, cl_norm, cl, cd])
            # print(frmstr.format(j, yle, chord, area, c_cl, ai, cl_norm, cl, cd))
        print(table)
    def print_strip_coefficients(self):
        from py2md.classes import MDTable
        from math import atan
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('Chord', '.5f')
        table.add_column('Area', '.5f')
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
                force += self.nffrc[pnl.lpid]
                momle += self.nfmom[pnl.lpid]
                rref = pnl.pnti-strp.pnti
                momle += rref**self.nffrc[pnl.lpid]
            cn = force*strp.nrmt/q/area
            ca = force.x/q/area
            cl = self.trlft[strp.lsid]/q/area
            cd = self.trdrg[strp.lsid]/q/area
            dw = -self.trwsh[strp.lsid]
            cmle = momle.y/q/area/chord
            rqc = Vector(-chord/4, 0.0, 0.0)
            momqc = momle+rqc**force
            cmqc = momqc.y/q/area/chord
            table.add_row([j, chord, area, cn, ca, cl, cd, dw, cmle, cmqc])
        print(table)
    def print_panel_forces(self):
        from py2md.classes import MDTable
        from math import atan, tan, radians
        table = MDTable()
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
            x = pnl.pnti.x
            y = pnl.pnti.y
            z = pnl.pnti.z
            area = pnl.area
            frc = self.nffrc[pnl.lpid]
            nfrc = frc*pnl.nrml
            cp = nfrc/area/q
            chord = pnl.crd
            alc = tan(radians(pnl.alpha))
            table.add_row([j, k, x, y, z, chord, alc, cp])
        print(table)
    def print_panel_near_field_results(self):
        from py2md.classes import MDTable
        from math import atan
        table = MDTable()
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
            gam = self.gam[j]
            Vx = self.nfvel[j].x
            Vy = self.nfvel[j].y
            Vz = self.nfvel[j].z
            lx = pnl.leni.x
            ly = pnl.leni.y
            lz = pnl.leni.z
            Fx = self.nffrc[j].x
            Fy = self.nffrc[j].y
            Fz = self.nffrc[j].z
            table.add_row([j, k, gam, Vx, Vy, Vz, lx, ly, lz, Fx, Fy, Fz])
        print(table)
    def __str__(self):
        from py2md.classes import MDTable
        from . import cfrm, dfrm, efrm
        outstr = '# Lattice Result '+self.name+'\n'
        table = MDTable()
        table.add_column('Alpha (deg)', 'g', data=[self.alpha])
        table.add_column('Beta (deg)', 'g', data=[self.beta])
        table.add_column('Speed', 'g', data=[self.speed])
        table.add_column('Rho', 'g', data=[self.rho])
        outstr += table._repr_markdown_()
        if self.gam is not None:
            table = MDTable()
            table.add_column('Cx', cfrm, data=[self.Cx])
            table.add_column('Cy', cfrm, data=[self.Cy])
            table.add_column('Cz', cfrm, data=[self.Cz])
            table.add_column('Cl', cfrm, data=[self.Cl])
            table.add_column('Cm', cfrm, data=[self.Cm])
            table.add_column('Cn', cfrm, data=[self.Cn])
            outstr += table._repr_markdown_()
            table = MDTable()
            table.add_column('CL', cfrm, data=[self.CL])
            table.add_column('CDi', dfrm, data=[self.CDi])
            table.add_column('CY', cfrm, data=[self.CY])
            outstr += table._repr_markdown_()
        if self.phi is not None:
            table = MDTable()
            table.add_column('CL_ff', cfrm, data=[self.CL_ff])
            table.add_column('CDi_ff', dfrm, data=[self.CDi_ff])
            table.add_column('CY_ff', cfrm, data=[self.CY_ff])
            table.add_column('e', efrm, data=[self.e])
            outstr += table._repr_markdown_()
        return outstr
    def __repr__(self):
        return f'<LatticeResult: {self.name}>'
    def _repr_markdown_(self):
        return self.__str__()

class PanelResult(object):
    pnl = None
    gam = None
    vel = None
    frc = None
    mom = None
    def __init__(self, pnl):
        self.pnl = pnl
    def set_gam(self, gam: float):
        self.gam = gam
    def set_velocity(self, vel: Vector):
        self.vel = vel
    def set_force(self, frc: Vector):
        self.frc = frc
    def set_moment(self, mom: Vector):
        self.mom = mom

class StripResult(object):
    strp = None
    phi = None
    nfwsh = None
    nflft = None
    nfdrg = None
    nfmom = None
    pnlres = None
    def __init__(self, strp):
        self.strp = strp
        self.pnlres = []
    def add_panel_result(self, pnlres):
        self.pnlres.append(pnlres)
