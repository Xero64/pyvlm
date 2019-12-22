from pygeom.geom3d import Point, Coordinate, Vector, jhat, ihat
from pygeom.matrix3d import zero_matrix_vector, elementwise_multiply
from pygeom.matrix3d import elementwise_cross_product, elementwise_dot_product
from numpy.matlib import zeros, matrix
from numpy import multiply
from math import radians, cos, sin
from matplotlib.pyplot import figure
from .latticesystem import LatticeSystem

class LatticeResult(object):
    name = None
    sys = None
    rho = None
    speed = None
    alpha = None
    beta = None
    pbo2V = None
    qco2V = None
    rbo2V = None
    ctrls = None
    rcg = None
    _acs = None
    _scs = None
    _wcs = None
    _vfs = None
    _qfs = None
    _ofs = None
    _gamma = None
    _avv = None
    _afv = None
    _arm = None
    _phi = None
    _bvv = None
    _brm = None
    _nfres = None
    _trres = None
    _pdres = None
    _stgam = None
    _stres = None
    _ctgamp = None
    _ctresp = None
    _ctgamn = None
    _ctresn = None
    def __init__(self, name: str, sys: LatticeSystem):
        self.name = name
        self.sys = sys
        self.initialise()
    def initialise(self):
        self.rho = 1.0
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
    def reset(self):
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None
    def set_density(self, rho: float=None):
        if rho is not None:
            self.rho = rho
        self.reset()
    def set_state(self, speed: float=None, alpha: float=None, beta: float=None,
                  pbo2V: float=None, qco2V: float=None, rbo2V: float=None):
        if speed is not None:
            self.speed = speed
        if alpha is not None:
            self.alpha = alpha
        if beta is not None:
            self.beta = beta
        # print(f'name = {self.name}')
        if pbo2V is not None:
            self.pbo2V = pbo2V
            # print(f'pbo2V = {self.pbo2V}')
        if qco2V is not None:
            self.qco2V = qco2V
            # print(f'qco2V = {self.qco2V}')
        if rbo2V is not None:
            self.rbo2V = rbo2V
            # print(f'rbo2V = {self.rbo2V}')
        self.reset()
    def set_controls(self, **kwargs):
        for control in kwargs:
            self.ctrls[control] = kwargs[control]
        self.reset()
    def set_cg(self, rcg: Point):
        self.rcg = rcg
        self.reset()
    @property
    def acs(self):
        if self._acs is None:
            pnt = self.sys.rref
            cosal, sinal = trig_angle(self.alpha)
            cosbt, sinbt = trig_angle(self.beta)
            dirx = Vector(cosbt*cosal, -sinbt, cosbt*sinal)
            diry = Vector(sinbt*cosal, cosbt, sinbt*sinal)
            dirz = Vector(-sinal, 0.0, cosal)
            self._acs = Coordinate(pnt, dirx, diry, dirz)
        return self._acs
    @property
    def scs(self):
        if self._scs is None:
            pnt = self.sys.rref
            cosal, sinal = trig_angle(self.alpha)
            dirx = Vector(cosal, 0.0, sinal)
            diry = Vector(0.0, 1.0, 0.0)
            dirz = Vector(-sinal, 0.0, cosal)
            self._scs = Coordinate(pnt, dirx, diry, dirz)
        return self._scs
    @property
    def wcs(self):
        if self._wcs is None:
            pnt = self.sys.rref
            dirx = -1.0*self.acs.dirx
            diry = self.acs.diry
            dirz = -1.0*self.acs.dirz
            self._wcs = Coordinate(pnt, dirx, diry, dirz)
        return self._wcs
    def set_gamma(self, gam: list):
        self._gamma = matrix([gam], dtype=float).transpose()
    def set_phi(self, phi: list):
        if len(phi) != len(self.sys.strps):
            raise Exception('The length of phi must equal the number of strips.')
        self._phi = matrix([phi], dtype=float).transpose()
    def set_lift_distribution(self, l: list, rho: float, speed: float):
        if len(l) != len(self.sys.strps):
            raise Exception('The length of l must equal the number of strips.')
        phi = [li/rho/speed for li in l]
        self.set_density(rho)
        self.set_state(speed)
        self.set_phi(phi)
    @property
    def gamma(self):
        if self._gamma is None:
            self._gamma = vector_matrix_dot(self.sys.gam[:, 0], self.vfs)
            self._gamma += vector_matrix_dot(self.sys.gam[:, 1], self.ofs)
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
                    self._gamma += vector_matrix_dot(self.sys.gam[:, indv], self.vfs*ctrlrad)
                    self._gamma += vector_matrix_dot(self.sys.gam[:, indo], self.ofs*ctrlrad)
        return self._gamma
    def galpha(self):
        vfs = dufsdal_sa(radians(self.alpha))*self.speed
        vfs, ofs = dalpha(self.speed, self.alpha, self.beta,
                          self.pbo2V, self.qco2V, self.rbo2V,
                          self.sys.bref, self.sys.cref)
        gmat = vector_matrix_dot(self.sys.gam[:, 0], vfs)
        gmat += vector_matrix_dot(self.sys.gam[:, 1], ofs)
        return gmat
    def gbeta(self):
        vfs, ofs = dbeta(self.speed, self.alpha, self.beta,
                         self.pbo2V, self.qco2V,
                         self.sys.bref, self.sys.cref)
        gmat = vector_matrix_dot(self.sys.gam[:, 0], vfs)
        gmat += vector_matrix_dot(self.sys.gam[:, 1], ofs)
        return gmat
    def gpb02V(self):
        ofs = dpbo2V(self.speed, self.alpha, self.beta, self.sys.bref)
        gmat = vector_matrix_dot(self.sys.gam[:, 1], ofs)
        return gmat
    def gqc02V(self):
        ofs = dqco2V(self.speed, self.alpha, self.beta, self.sys.cref)
        gmat = vector_matrix_dot(self.sys.gam[:, 1], ofs)
        return gmat
    def grb02V(self):
        ofs = drbo2V(self.speed, self.alpha, self.sys.bref)
        gmat = vector_matrix_dot(self.sys.gam[:, 1], ofs)
        return gmat
    def gctrlp_single(self, control: str):
        indv = self.sys.ctrls[control][0]
        gmat = vector_matrix_dot(self.sys.gam[:, indv], self.vfs)
        indo = self.sys.ctrls[control][1]
        gmat += vector_matrix_dot(self.sys.gam[:, indo], self.ofs)
        return gmat
    def gctrlp(self, control: str=''):
        gmats = {}
        for control in self.sys.ctrls:
            gmats[control] = self.gctrlp_single(control)
        return gmats
    def gctrln_single(self, control: str):
        indv = self.sys.ctrls[control][2]
        gmat = vector_matrix_dot(self.sys.gam[:, indv], self.vfs)
        indo = self.sys.ctrls[control][3]
        gmat += vector_matrix_dot(self.sys.gam[:, indo], self.ofs)
        return gmat
    def gctrln(self):
        gmats = {}
        for control in self.sys.ctrls:
            gmats[control] = self.gctrln_single(control)
        return gmats 
    @property
    def phi(self):
        if self._phi is None:
            num = len(self.sys.strps)
            self._phi = zeros((num, 1))
            for strp in self.sys.strps:
                i = strp.lsid
                for pnl in strp.pnls:
                    self._phi[i, 0] += self.gamma[pnl.lpid, 0]
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
            self._vfs = self.acs.dirx*self.speed
        return self._vfs
    def calc_ofs(self, pbo2V: float, qco2V: float, rbo2V: float):
        p = pbo2V*2*self.speed/self.sys.bref
        q = qco2V*2*self.speed/self.sys.cref
        r = rbo2V*2*self.speed/self.sys.bref
        rotvl = Vector(p, q, r)
        return self.wcs.vector_to_global(rotvl)
    @property
    def ofs(self):
        if self._ofs is None:
            self._ofs = self.calc_ofs(self.pbo2V, self.qco2V, self.rbo2V)
        return self._ofs
    @property
    def qfs(self):
        if self._qfs is None:
            self._qfs = self.rho*self.speed**2/2
        return self._qfs
    @property
    def avv(self):
        if self._avv is None:
            num = len(self.sys.pnls)
            self._avv = zero_matrix_vector((num, 1))
            for pnl in self.sys.pnls:
                i = pnl.lpid
                if pnl.noload:
                    self._avv[i, 0] = Vector(0.0, 0.0, 0.0)
                else:
                    self._avv[i, 0] = self.vfs-self.ofs**self.arm[i, 0]
        return self._avv
    @property
    def bvv(self):
        if self._bvv is None:
            num = len(self.sys.strps)
            self._bvv = zero_matrix_vector((num, 1))
            for strp in self.sys.strps:
                i = strp.lsid
                if strp.noload:
                    self._bvv[i, 0] = Vector(0.0, 0.0, 0.0)
                else:
                    self._bvv[i, 0] = self.vfs-self.ofs**self.brm[i, 0]
        return self._bvv
    @property
    def arm(self):
        if self._arm is None:
            num = len(self.sys.pnls)
            self._arm = zero_matrix_vector((num, 1))
            for pnl in self.sys.pnls:
                i = pnl.lpid
                if pnl.noload:
                    self._arm[i, 0] = Vector(0.0, 0.0, 0.0)
                else:
                    self._arm[i, 0] = pnl.pnti-self.rcg
        return self._arm
    @property
    def brm(self):
        if self._brm is None:
            num = len(self.sys.strps)
            self._brm = zero_matrix_vector((num, 1))
            for strp in self.sys.strps:
                i = strp.lsid
                if strp.noload:
                    self._brm[i, 0] = Vector(0.0, 0.0, 0.0)
                else:
                    self._brm[i, 0] = strp.pntq-self.rcg
        return self._brm
    @property
    def afv(self):
        if self._afv is None:
            num = len(self.sys.pnls)
            self._afv = zero_matrix_vector((num, 1))
            for pnl in self.sys.pnls:
                i = pnl.lpid
                self._afv[i, 0] = pnl.induced_force(self.avv[i, 0])
        return self._afv
    @property
    def nfres(self):
        if self._nfres is None:
            self._nfres = GammaResult(self, self.gamma)
        return self._nfres
    @property
    def trres(self):
        if self._trres is None:
            self._trres = PhiResult(self, self.phi)
        return self._trres
    @property
    def pdres(self):
        if self._pdres is None:
            self._pdres = ParasiticDragResult(self)
        return self._pdres
    @property
    def stgam(self):
        if self._stgam is None:
            self._stgam = {}
            self._stgam['alpha'] = self.galpha()
            self._stgam['beta'] = self.gbeta()
            self._stgam['pbo2V'] = self.gpb02V()
            self._stgam['qco2V'] = self.gqc02V()
            self._stgam['rbo2V'] = self.grb02V()
        return self._stgam
    @property
    def stres(self):
        if self._stres is None:
            self._stres = {}
            self._stres['alpha'] = GammaResult(self, self.stgam['alpha'])
            self._stres['beta'] = GammaResult(self, self.stgam['beta'])
            self._stres['pbo2V'] = GammaResult(self, self.stgam['pbo2V'])
            self._stres['qco2V'] = GammaResult(self, self.stgam['qco2V'])
            self._stres['rbo2V'] = GammaResult(self, self.stgam['rbo2V'])
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
    def plot_panel_near_field_velocities(self, ax=None, component=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        py = [pnl.pnti.y for pnl in self.sys.pnls]
        if component is None or component == 'x':
            vx = [vel.x for vel in self.nfres.nfvel.transpose().tolist()[0]]
            ax.plot(py, vx, label=self.name+' Velocity X')
        if component is None or component == 'y':
            vy = [vel.y for vel in self.nfres.nfvel.transpose().tolist()[0]]
            ax.plot(py, vy, label=self.name+' Velocity Y')
        if component is None or component == 'z':
            vz = [vel.z for vel in self.nfres.nfvel.transpose().tolist()[0]]
            ax.plot(py, vz, label=self.name+' Velocity Z')
        ax.legend()
        return ax
    def plot_phi_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.sys.srfcs[0].strpy, self.phi, label=self.name)
        ax.legend()
        return ax
    def plot_trefftz_lift_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        for srfc in self.sys.srfcs:
            y = []
            l = []
            for strp in srfc.strps:
                lsid = strp.lsid
                if strp.dyt != 0.0:
                    y.append(strp.pnti.y)
                    l.append(self.trres.trfrc.z[lsid, 0]/strp.dyt)
            if len(l) > 0:
                label = self.name+' for '+srfc.name
                ax.plot(y, l, label=label)
        ax.legend()
        return ax
    def plot_trefftz_yforce_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        for srfc in self.sys.srfcs:
            z = []
            f = []
            for strp in srfc.strps:
                lsid = strp.lsid
                if strp.dzt != 0.0:
                    z.append(strp.pnti.z)
                    f.append(self.trres.trfrc.y[lsid, 0]/strp.dzt)
            if len(f) > 0:
                label = self.name+' for '+srfc.name
                ax.plot(f, z, label=label)
        ax.legend()
        return ax
    def plot_trefftz_drag_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        for srfc in self.sys.srfcs:
            y = []
            d = []
            for strp in srfc.strps:
                lsid = strp.lsid
                y.append(strp.pnti.y)
                d.append(self.trres.trfrc.x[lsid, 0]/strp.dst)
            if len(d) > 0:
                label = self.name+' for '+srfc.name
                ax.plot(y, d, label=label)
        ax.legend()
        return ax
    def plot_trefftz_wash_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        for srfc in self.sys.srfcs:
            y = []
            w = []
            for strp in srfc.strps:
                lsid = strp.lsid
                y.append(strp.pnti.y)
                w.append(self.trres.trwsh[lsid, 0])
            if len(w) > 0:
                label = self.name+' for '+srfc.name
                ax.plot(y, w, label=label)
        ax.legend()
        return ax
    def to_result(self, name: str=''):
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
    def surface_loads(self):
        from py2md.classes import MDTable
        from math import atan
        table = MDTable()
        table.add_column('Name', 's')
        table.add_column('Fx', '.3f')
        table.add_column('Fy', '.3f')
        table.add_column('Fz', '.3f')
        table.add_column('Mx', '.3f')
        table.add_column('My', '.3f')
        table.add_column('Mz', '.3f')
        for srfc in self.sys.srfcs:
            ind = srfc.pnli
            frc = self.nfres.nffrc[ind, 0].sumall()
            mom = self.nfres.nfmom[ind, 0].sumall()
            table.add_row([srfc.name, frc.x, frc.y, frc.z, mom.x, mom.y, mom.z])
        frc = self.nfres.nffrc.sumall()
        mom = self.nfres.nfmom.sumall()
        table.add_row(['Total', frc.x, frc.y, frc.z, mom.x, mom.y, mom.z])
        return table
    @property
    def strip_forces(self):
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
                nrmfrc += self.nfres.nffrc[pnl.lpid, 0]/pnl.area*pnl.crd
            c_cl = nrmfrc.z/q
            ai = -self.trres.trwsh[strp.lsid, 0]/self.speed
            cd = self.trres.trfrc.x[strp.lsid, 0]/q/area
            cy = nrmfrc*self.acs.diry/q/chord
            cz = nrmfrc*self.acs.dirz/q/chord
            cf = Vector(0.0, cy, cz)
            cl_norm = cf*strp.nrmt
            cl = cl_norm
            table.add_row([j, yle, chord, area, c_cl, ai, cl_norm, cl, cd])
        return table
    @property
    def strip_coefficients(self):
        from py2md.classes import MDTable
        from math import atan
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('Chord', '.4f')
        table.add_column('Area', '.6f')
        table.add_column('cn', '.5f') # Strip Normal Load
        table.add_column('ca', '.5f') # String Axial Load
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
                force += self.nfres.nffrc[pnl.lpid, 0]
                rref = pnl.pnti-strp.pnti
                momle += rref**self.nfres.nffrc[pnl.lpid, 0]
            cn = force*strp.nrmt/q/area
            ca = force.x/q/area
            cl = force*self.acs.dirz/q/area
            cd = force*self.acs.dirx/q/area
            dw = -self.trres.trwsh[strp.lsid, 0]/self.speed
            cmle = momle.y/q/area/chord
            rqc = Vector(-chord/4, 0.0, 0.0)
            momqc = momle+rqc**force
            cmqc = momqc.y/q/area/chord
            table.add_row([j, chord, area, cn, ca, cl, cd, dw, cmle, cmqc])
        return table
    @property
    def panel_forces(self):
        from py2md.classes import MDTable
        from math import atan, tan
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
            x = pnl.pntg.x
            y = pnl.pntg.y
            z = pnl.pntg.z
            area = pnl.area
            frc = self.nfres.nffrc[pnl.lpid, 0]
            nfrc = frc*pnl.nrml
            cp = nfrc/area/q
            chord = pnl.crd
            alc = tan(radians(pnl.alpha))
            table.add_row([j, k, x, y, z, chord, alc, cp])
        return table
    @property
    def panel_near_field_results(self):
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
            gamma = self.gamma[j, 0]
            Vx = self.nfres.nfvel[j, 0].x
            Vy = self.nfres.nfvel[j, 0].y
            Vz = self.nfres.nfvel[j, 0].z
            lx = pnl.leni.x
            ly = pnl.leni.y
            lz = pnl.leni.z
            Fx = self.nfres.nffrc[j, 0].x
            Fy = self.nfres.nffrc[j, 0].y
            Fz = self.nfres.nffrc[j, 0].z
            table.add_row([j, k, gamma, Vx, Vy, Vz, lx, ly, lz, Fx, Fy, Fz])
        return table
    @property
    def stability_derivatives(self):
        from py2md.classes import MDTable, MDHeading, MDReport
        from . import sfrm
        report = MDReport()
        heading = MDHeading('Stability Derivatives', 1)
        report.add_object(heading)
        table = MDTable()
        table.add_column('CLa', sfrm, data=[self.stres['alpha'].CL])
        table.add_column('CYa', sfrm, data=[self.stres['alpha'].CY])
        table.add_column('Cla', sfrm, data=[self.stres['alpha'].Cl])
        table.add_column('Cma', sfrm, data=[self.stres['alpha'].Cm])
        table.add_column('Cna', sfrm, data=[self.stres['alpha'].Cn])
        report.add_object(table)
        table = MDTable()
        table.add_column('CLb', sfrm, data=[self.stres['beta'].CL])
        table.add_column('CYb', sfrm, data=[self.stres['beta'].CY])
        table.add_column('Clb', sfrm, data=[self.stres['beta'].Cl])
        table.add_column('Cmb', sfrm, data=[self.stres['beta'].Cm])
        table.add_column('Cnb', sfrm, data=[self.stres['beta'].Cn])
        report.add_object(table)
        table = MDTable()
        table.add_column('CLp', sfrm, data=[self.stres['pbo2V'].CL])
        table.add_column('CYp', sfrm, data=[self.stres['pbo2V'].CY])
        table.add_column('Clp', sfrm, data=[self.stres['pbo2V'].Cl])
        table.add_column('Cmp', sfrm, data=[self.stres['pbo2V'].Cm])
        table.add_column('Cnp', sfrm, data=[self.stres['pbo2V'].Cn])
        report.add_object(table)
        table = MDTable()
        table.add_column('CLq', sfrm, data=[self.stres['qco2V'].CL])
        table.add_column('CYq', sfrm, data=[self.stres['qco2V'].CY])
        table.add_column('Clq', sfrm, data=[self.stres['qco2V'].Cl])
        table.add_column('Cmq', sfrm, data=[self.stres['qco2V'].Cm])
        table.add_column('Cnq', sfrm, data=[self.stres['qco2V'].Cn])
        report.add_object(table)
        table = MDTable()
        table.add_column('CLr', sfrm, data=[self.stres['rbo2V'].CL])
        table.add_column('CYr', sfrm, data=[self.stres['rbo2V'].CY])
        table.add_column('Clr', sfrm, data=[self.stres['rbo2V'].Cl])
        table.add_column('Cmr', sfrm, data=[self.stres['rbo2V'].Cm])
        table.add_column('Cnr', sfrm, data=[self.stres['rbo2V'].Cn])
        report.add_object(table)
        return report
    @property
    def control_derivatives(self):
        from py2md.classes import MDTable, MDHeading, MDReport
        from . import sfrm
        report = MDReport()
        heading = MDHeading('Control Derivatives', 1)
        report.add_object(heading)
        for control in self.ctrls:
            letter = control[0]
            heading = MDHeading(f'{control.capitalize()} Derivatives', 2)
            report.add_object(heading)
            table = MDTable()
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
            report.add_object(table)
        return report
    def __str__(self):
        from py2md.classes import MDTable
        from . import cfrm, dfrm, efrm
        outstr = '# Lattice Result '+self.name+' for '+self.sys.name+'\n'
        table = MDTable()
        table.add_column('Alpha (deg)', cfrm, data=[self.alpha])
        table.add_column('Beta (deg)', cfrm, data=[self.beta])
        table.add_column('Speed', cfrm, data=[self.speed])
        table.add_column('Rho', cfrm, data=[self.rho])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('pb/2V (rad)', cfrm, data=[self.pbo2V])
        table.add_column('qc/2V (rad)', cfrm, data=[self.qco2V])
        table.add_column('rb/2V (rad)', cfrm, data=[self.rbo2V])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('xcg', '.5f', data=[self.rcg.x])
        table.add_column('ycg', '.5f', data=[self.rcg.y])
        table.add_column('zcg', '.5f', data=[self.rcg.z])
        outstr += table._repr_markdown_()
        if len(self.ctrls) > 0:
            table = MDTable()
            for control in self.ctrls:
                ctrl = self.ctrls[control]
                control = control.capitalize()
                table.add_column(f'{control} (deg)', cfrm, data=[ctrl])
            outstr += table._repr_markdown_()
        if self.sys.cdo != 0.0:
            table = MDTable()
            table.add_column('CDo', dfrm, data=[self.pdres.CDo])
            table.add_column('CYo', cfrm, data=[self.pdres.CY])
            table.add_column('CLo', cfrm, data=[self.pdres.CL])
            table.add_column('Clo', cfrm, data=[self.pdres.Cl])
            table.add_column('Cmo', cfrm, data=[self.pdres.Cm])
            table.add_column('Cno', cfrm, data=[self.pdres.Cn])
            outstr += table._repr_markdown_()
        if self.gamma is not None:
            table = MDTable()
            table.add_column('Cx', cfrm, data=[self.nfres.Cx])
            table.add_column('Cy', cfrm, data=[self.nfres.Cy])
            table.add_column('Cz', cfrm, data=[self.nfres.Cz])
            outstr += table._repr_markdown_()
            table = MDTable()
            table.add_column('CDi', dfrm, data=[self.nfres.CDi])
            table.add_column('CY', cfrm, data=[self.nfres.CY])
            table.add_column('CL', cfrm, data=[self.nfres.CL])
            table.add_column('Cl', cfrm, data=[self.nfres.Cl])
            table.add_column('Cm', cfrm, data=[self.nfres.Cm])
            table.add_column('Cn', cfrm, data=[self.nfres.Cn])
            table.add_column('e', efrm, data=[self.nfres.e])
            if self.sys.cdo != 0.0:
                lod = self.nfres.CL/(self.pdres.CDo+self.nfres.CDi)
                table.add_column('L/D', '.5g', data=[lod])
            outstr += table._repr_markdown_()
        if self.phi is not None:
            table = MDTable()
            table.add_column('CDi_ff', dfrm, data=[self.trres.CDi])
            table.add_column('CY_ff', cfrm, data=[self.trres.CY])
            table.add_column('CL_ff', cfrm, data=[self.trres.CL])
            # table.add_column('Cl_ff', cfrm, data=[self.trres.Cl])
            # table.add_column('Cm_ff', cfrm, data=[self.trres.Cm])
            # table.add_column('Cn_ff', cfrm, data=[self.trres.Cn])
            table.add_column('e', efrm, data=[self.trres.e])
            if self.sys.cdo != 0.0:
                lod_ff = self.trres.CL/(self.pdres.CDo+self.trres.CDi)
                table.add_column('L/D_ff', '.5g', data=[lod_ff])
            outstr += table._repr_markdown_()
        return outstr
    def __repr__(self):
        return f'<LatticeResult: {self.name}>'
    def _repr_markdown_(self):
        return self.__str__()

def trig_angle(angle: float):
    '''Calculates cos(angle) and sin(angle) with angle in degrees.'''
    angrad = radians(angle)
    cosang = cos(angrad)
    sinang = sin(angrad)
    return cosang, sinang

def dufsdal(cosal: float, sinal: float, cosbt: float, sinbt: float):
    return Vector(-sinal*cosbt, 0.0, cosal*cosbt)

def dufsdbt(cosal: float, sinal: float, cosbt: float, sinbt: float):
    return Vector(-sinbt*cosal, -cosbt, -sinal*sinbt)

def dufsdal_sa(alpha: float):
    return Vector(-alpha, 0.0, 1.0)

def dufsdbt_sa(beta: float):
    return Vector(-beta, -1.0, 0.0)

def dalpha(V: float, al: float, bt: float, p: float, q: float, r: float, b: float, c: float):
    cosal, sinal = trig_angle(al)
    cosbt, sinbt = trig_angle(bt)
    val = Vector(-V*sinal*cosbt, 0.0, V*cosal*cosbt)
    oal = Vector(2*V*(q*sinal*sinbt/c + p*sinal*cosbt/b + r*cosal/b), 0.0,
                 -2*V*(q*sinbt*cosal/c + p*cosal*cosbt/b - r*sinal/b))
    return val, oal

def dbeta(V: float, al: float, bt: float, p: float, q: float, b: float, c: float):
    cosal, sinal = trig_angle(al)
    cosbt, sinbt = trig_angle(bt)
    vbt = Vector(-V*sinbt*cosal, -V*cosbt, -V*sinal*sinbt)
    obt = Vector(2*V*(-q*cosbt/c + p*sinbt/b)*cosal,
                 -2*V*(q*sinbt/c + p*cosbt/b),
                 -2*V*(q*cosbt/c - p*sinbt/b)*sinal)
    return vbt, obt

def dpbo2V(V: float, al: float, bt: float, b: float):
    cosal, sinal = trig_angle(al)
    cosbt, sinbt = trig_angle(bt)
    return Vector(-2*V*cosal*cosbt/b, -2*V*sinbt/b, -2*V*sinal*cosbt/b)

def dqco2V(V: float, al: float, bt: float, c: float):
    cosal, sinal = trig_angle(al)
    cosbt, sinbt = trig_angle(bt)
    return Vector(-2*V*sinbt*cosal/c, 2*V*cosbt/c, -2*V*sinal*sinbt/c)

def drbo2V(V: float, al: float, b: float):
    cosal, sinal = trig_angle(al)
    return Vector(2*V*sinal/b, 0.0, -2*V*cosal/b)

def vector_matrix_dot(mat: matrix, vec: Vector):
    outmat = zeros(mat.shape)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            outmat[i, j] += mat[i, j]*vec
    return outmat

class GammaResult(object):
    res = None
    gamma = None
    _rhogamma = None
    _nfvel = None
    _nffrc = None
    _nfmom = None
    _nffrctot = None
    _nfmomtot = None
    _Cx = None
    _Cy = None
    _Cz = None
    _CDi = None
    _CY = None
    _CL = None
    _e = None
    _Cl = None
    _Cm = None
    _Cn = None
    def __init__(self, res: LatticeResult, gamma: matrix):
        self.res = res
        self.gamma = gamma
    @property
    def rhogamma(self):
        if self._rhogamma is None:
            self._rhogamma = self.res.rho*self.gamma
        return self._rhogamma
    @property
    def nfvel(self):
        if self._nfvel is None:
            self._nfvel = self.res.sys.avg*self.gamma+self.res.avv
        return self._nfvel
    @property
    def nffrc(self):
        if self._nffrc is None:
            tmp = self.res.sys.afg*self.gamma+self.res.afv
            self._nffrc = elementwise_multiply(self.rhogamma, tmp)
        return self._nffrc
    @property
    def nfmom(self):
        if self._nfmom is None:
            self._nfmom = elementwise_cross_product(self.res.arm, self.nffrc)
        return self._nfmom
    @property
    def nffrctot(self):
        if self._nffrctot is None:
            self._nffrctot = self.nffrc.sumall()
        return self._nffrctot
    @property
    def nfmomtot(self):
        if self._nfmomtot is None:
            self._nfmomtot = self.nfmom.sumall()
        return self._nfmomtot
    @property
    def Cx(self):
        if self._Cx is None:
            self._Cx = self.nffrctot.x/self.res.qfs/self.res.sys.sref
            self._Cx = fix_zero(self._Cx)
        return self._Cx
    @property
    def Cy(self):
        if self._Cy is None:
            self._Cy = self.nffrctot.y/self.res.qfs/self.res.sys.sref
            self._Cy = fix_zero(self._Cy)
        return self._Cy
    @property
    def Cz(self):
        if self._Cz is None:
            self._Cz = self.nffrctot.z/self.res.qfs/self.res.sys.sref
            self._Cz = fix_zero(self._Cz)
        return self._Cz
    @property
    def CDi(self):
        if self._CDi is None:
            Di = self.res.acs.dirx*self.nffrctot
            self._CDi = Di/self.res.qfs/self.res.sys.sref
            self._CDi = fix_zero(self._CDi)
        return self._CDi
    @property
    def CY(self):
        if self._CY is None:
            Y = self.res.acs.diry*self.nffrctot
            self._CY = Y/self.res.qfs/self.res.sys.sref
            self._CY = fix_zero(self._CY)
        return self._CY
    @property
    def CL(self):
        if self._CL is None:
            L = self.res.acs.dirz*self.nffrctot
            self._CL = L/self.res.qfs/self.res.sys.sref
            self._CL = fix_zero(self._CL)
        return self._CL
    @property
    def Cl(self):
        if self._Cl is None:
            l = self.res.wcs.dirx*self.nfmomtot
            self._Cl = l/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cl = fix_zero(self._Cl)
        return self._Cl
    @property
    def Cm(self):
        if self._Cm is None:
            m = self.res.wcs.diry*self.nfmomtot
            self._Cm = m/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cm = fix_zero(self._Cm)
        return self._Cm
    @property
    def Cn(self):
        if self._Cn is None:
            n = self.res.wcs.dirz*self.nfmomtot
            self._Cn = n/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cn = fix_zero(self._Cn)
        return self._Cn
    @property
    def e(self):
        if self._e is None:
            if self.CDi <= 0.0:
                self._e = float('nan')
            elif self.CL == 0.0 and self.CY == 0.0:
                self._e = 0.0
            else:
                from math import pi
                self._e = (self.CL**2+self.CY**2)/pi/self.res.sys.ar/self.CDi
                self._e = fix_zero(self._e)
        return self._e

class PhiResult(object):
    res = None
    phi = None
    _trwsh = None
    _trfrc = None
    _trmom = None
    _trfrctot = None
    _trmomtot = None
    _CDi = None
    _CY = None
    _CL = None
    _Cl = None
    _Cm = None
    _Cn = None
    _e = None
    _lod = None
    def __init__(self, res: LatticeResult, phi: matrix):
        self.res = res
        self.phi = phi
    @property
    def trwsh(self):
        if self._trwsh is None:
            self._trwsh = self.res.sys.bvg*self.phi
        return self._trwsh
    @property
    def trfrc(self):
        if self._trfrc is None:
            from pygeom.matrix3d import MatrixVector
            x = self.res.rho*multiply(self.phi, self.res.sys.bdg*self.phi)
            y = self.res.rho*self.res.speed*multiply(self.phi, self.res.sys.byg)
            z = self.res.rho*self.res.speed*multiply(self.phi, self.res.sys.blg)
            self._trfrc = MatrixVector(x, y, z)
        return self._trfrc
    @property
    def trmom(self):
        if self._trmom is None:
            self._trmom = elementwise_cross_product(self.res.brm, self.trfrc)
        return self._trmom
    @property
    def trfrctot(self):
        if self._trfrctot is None:
            self._trfrctot = self.trfrc.sumall()
        return self._trfrctot
    @property
    def trmomtot(self):
        if self._trmomtot is None:
            self._trmomtot = self.trmom.sumall()
        return self._trmomtot
    @property
    def CDi(self):
        if self._CDi is None:
            Di = self.trfrctot.x
            self._CDi = Di/self.res.qfs/self.res.sys.sref
            self._CDi = fix_zero(self._CDi)
        return self._CDi
    @property
    def CY(self):
        if self._CY is None:
            Y = self.trfrctot.y
            self._CY = Y/self.res.qfs/self.res.sys.sref
            self._CY = fix_zero(self._CY)
        return self._CY
    @property
    def CL(self):
        if self._CL is None:
            L = self.trfrctot.z
            self._CL = L/self.res.qfs/self.res.sys.sref
            self._CL = fix_zero(self._CL)
        return self._CL
    @property
    def Cl(self):
        if self._Cl is None:
            l = -self.trmomtot.x
            self._Cl = l/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cl = fix_zero(self._Cl)
        return self._Cl
    @property
    def Cm(self):
        if self._Cm is None:
            m = self.trmomtot.y
            self._Cm = m/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cm = fix_zero(self._Cm)
        return self._Cm
    @property
    def Cn(self):
        if self._Cn is None:
            n = -self._trmomtot.z
            self._Cn = n/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cn = fix_zero(self._Cn)
        return self._Cn
    @property
    def e(self):
        if self._e is None:
            if self.CDi == 0.0:
                if self.CL == 0.0 and self.CY == 0.0:
                    self._e = 0.0
                else:
                    self._e = float('nan')
            else:
                from math import pi
                self._e = (self.CL**2+self.CY**2)/pi/self.res.sys.ar/self.CDi
                self._e = fix_zero(self._e)
        return self._e

class ParasiticDragResult(object):
    res = None
    _pdfrc = None
    _pdmom = None
    _pdfrctot = None
    _pdmomtot = None
    _CDo = None
    _CY = None
    _CL = None
    _Cl = None
    _Cm = None
    _Cn = None
    def __init__(self, res: LatticeResult):
        self.res = res
    @property
    def pdfrc(self):
        if self._pdfrc is None:
            dynpr = (self.res.rho/2)*elementwise_dot_product(self.res.bvv, self.res.bvv)
            self._pdfrc = elementwise_multiply(self.res.sys.bda, dynpr*self.res.acs.dirx)
        return self._pdfrc
    @property
    def pdmom(self):
        if self._pdmom is None:
            self._pdmom = elementwise_cross_product(self.res.brm, self.pdfrc)
        return self._pdmom
    @property
    def pdfrctot(self):
        if self._pdfrctot is None:
            self._pdfrctot = self.pdfrc.sumall()
        return self._pdfrctot
    @property
    def pdmomtot(self):
        if self._pdmomtot is None:
            self._pdmomtot = self.pdmom.sumall()
        return self._pdmomtot
    @property
    def CDo(self):
        if self._CDo is None:
            Do = self.res.acs.dirx*self.pdfrctot
            self._CDo = Do/self.res.qfs/self.res.sys.sref
            self._CDo = fix_zero(self._CDo)
        return self._CDo
    @property
    def CY(self):
        if self._CY is None:
            Y = self.res.acs.diry*self.pdfrctot
            self._CY = Y/self.res.qfs/self.res.sys.sref
            self._CY = fix_zero(self._CY)
        return self._CY
    @property
    def CL(self):
        if self._CL is None:
            L = self.res.acs.dirz*self.pdfrctot
            self._CL = L/self.res.qfs/self.res.sys.sref
            self._CL = fix_zero(self._CL)
        return self._CL
    @property
    def Cl(self):
        if self._Cl is None:
            l = self.res.wcs.dirx*self.pdmomtot
            self._Cl = l/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cl = fix_zero(self._Cl)
        return self._Cl
    @property
    def Cm(self):
        if self._Cm is None:
            m = self.res.wcs.diry*self.pdmomtot
            self._Cm = m/self.res.qfs/self.res.sys.sref/self.res.sys.cref
            self._Cm = fix_zero(self._Cm)
        return self._Cm
    @property
    def Cn(self):
        if self._Cn is None:
            n = self.res.wcs.dirz*self.pdmomtot
            self._Cn = n/self.res.qfs/self.res.sys.sref/self.res.sys.bref
            self._Cn = fix_zero(self._Cn)
        return self._Cn

def fix_zero(value: float, tol: float=1e-12):
    if abs(value) < tol:
        value = 0.0
    return value

def latticeresult_from_json(lsys: LatticeSystem, resdata: dict):
    name = resdata['name']
    if 'inherit' in resdata:
        inherit = resdata['inherit']
        if inherit in lsys.results:
            lres = lsys.results[inherit].to_result(name=name)
    else:
        lres = LatticeResult(name, lsys)
    for key in resdata:
        if key == 'name':
            continue
        elif key == 'inherit':
            continue
        elif key == 'density':
            rho = resdata['density']
            lres.set_density(rho=rho)
        elif key == 'speed':
            speed = resdata['speed']
            lres.set_state(speed=speed)
        elif key ==  'alpha':
            alpha = resdata['alpha']
            lres.set_state(alpha=alpha)
        elif key ==  'beta':
            beta = resdata['beta']
            lres.set_state(beta=beta)
        elif key ==  'pbo2V':
            pbo2V = resdata['pbo2V']
            lres.set_state(pbo2V=pbo2V)
        elif key ==  'qco2V':
            qco2V = resdata['qco2V']
            lres.set_state(qco2V=qco2V)
        elif key ==  'rbo2V':
            rbo2V = resdata['rbo2V']
            lres.set_state(rbo2V=rbo2V)
        elif key in lres.ctrls:
            lres.ctrls[key] = resdata[key]
        elif key == 'rcg':
            rcgdata = resdata[key]
            rcg = Point(rcgdata['x'], rcgdata['y'], rcgdata['z'])
            lres.set_cg(rcg)
    lsys.results[name] = lres
    return lres