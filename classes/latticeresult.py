from pygeom.geom3d import Coordinate, Vector, jhat
from pygeom.matrixgeom3d import zero_matrix_vector, elementwise_multiply
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
    _acs = None
    _scs = None
    _wcs = None
    _vfs = None
    _qfs = None
    _ofs = None
    _gamma = None
    _avv = None
    _afv = None
    _amv = None
    _phi = None
    _nfres = None
    _trres = None
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
        if pbo2V is not None:
            self.pbo2V = pbo2V
        if qco2V is not None:
            self.qco2V = qco2V
        if rbo2V is not None:
            self.rbo2V = rbo2V
        self.reset()
    def set_controls(self, **kwargs):
        for control in kwargs:
            self.ctrls[control] = kwargs[control]
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
        self._gamma = matrix(gam, dtype=float).transpose()
    def set_phi(self, phi: list):
        if len(phi) != len(self.sys.strps):
            raise Exception('The length of phi must equal the number of strips.')
        self._phi = matrix(phi, dtype=float).transpose()
    def set_lift_distribution(self, l: list, rho: float, speed: float):
        if len(l) != len(self.sys.strps):
            raise Exception('The length of l must equal the number of strips.')
        phi = [li/rho/speed for li in l]
        self.set_phi(phi)
        self.set_density(rho)
        self.set_state(speed)
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
            for i, strp in enumerate(self.sys.strps):
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
            for i, pnl in enumerate(self.sys.pnls):
                if pnl.noload:
                    self._avv[i, 0] = Vector(0.0, 0.0, 0.0)
                else:
                    rpi = self.sys.pnls[i].pnti-self.sys.rref
                    self._avv[i, 0] = self.vfs+self.ofs**rpi
        return self._avv
    @property
    def afv(self):
        if self._afv is None:
            num = len(self.sys.pnls)
            self._afv = zero_matrix_vector((num, 1))
            for i, pnl in enumerate(self.sys.pnls):
                self._afv[i, 0] = pnl.induced_force(self.avv[i, 0])
        return self._afv
    @property
    def amv(self):
        if self._amv is None:
            num = len(self.sys.pnls)
            self._amv = zero_matrix_vector((num, 1))
            for i, pnl in enumerate(self.sys.pnls):
                self._amv[i, 0] = (pnl.pnti-self.sys.rref)**self.afv[i, 0]
        return self._amv
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
                    l.append(self.trres.trlft[lsid, 0]/strp.dyt)
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
                    f.append(self.trres.trfrc[lsid, 0]/strp.dzt)
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
                d.append(self.trres.trdrg[lsid, 0]/strp.dst)
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
            cd = self.trres.trdrg[strp.lsid, 0]/q/area
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
        return table._repr_markdown_()
    @property
    def stability_derivatives(self):
        from py2md.classes import MDTable
        from . import sfrm
        outstr = '\n# Stability Derivatives\n'
        table = MDTable()
        table.add_column('CLa', sfrm, data=[self.stres['alpha'].CL])
        table.add_column('CYa', sfrm, data=[self.stres['alpha'].CY])
        table.add_column('Cla', sfrm, data=[self.stres['alpha'].Cl])
        table.add_column('Cma', sfrm, data=[self.stres['alpha'].Cm])
        table.add_column('Cna', sfrm, data=[self.stres['alpha'].Cn])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('CLb', sfrm, data=[self.stres['beta'].CL])
        table.add_column('CYb', sfrm, data=[self.stres['beta'].CY])
        table.add_column('Clb', sfrm, data=[self.stres['beta'].Cl])
        table.add_column('Cmb', sfrm, data=[self.stres['beta'].Cm])
        table.add_column('Cnb', sfrm, data=[self.stres['beta'].Cn])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('CLp', sfrm, data=[self.stres['pbo2V'].CL])
        table.add_column('CYp', sfrm, data=[self.stres['pbo2V'].CY])
        table.add_column('Clp', sfrm, data=[self.stres['pbo2V'].Cl])
        table.add_column('Cmp', sfrm, data=[self.stres['pbo2V'].Cm])
        table.add_column('Cnp', sfrm, data=[self.stres['pbo2V'].Cn])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('CLq', sfrm, data=[self.stres['qco2V'].CL])
        table.add_column('CYq', sfrm, data=[self.stres['qco2V'].CY])
        table.add_column('Clq', sfrm, data=[self.stres['qco2V'].Cl])
        table.add_column('Cmq', sfrm, data=[self.stres['qco2V'].Cm])
        table.add_column('Cnq', sfrm, data=[self.stres['qco2V'].Cn])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('CLr', sfrm, data=[self.stres['rbo2V'].CL])
        table.add_column('CYr', sfrm, data=[self.stres['rbo2V'].CY])
        table.add_column('Clr', sfrm, data=[self.stres['rbo2V'].Cl])
        table.add_column('Cmr', sfrm, data=[self.stres['rbo2V'].Cm])
        table.add_column('Cnr', sfrm, data=[self.stres['rbo2V'].Cn])
        outstr += table._repr_markdown_()
        return outstr
    @property
    def control_derivatives(self):
        from py2md.classes import MDTable
        from . import sfrm
        outstr = '\n# Control Derivatives\n'
        for control in self.ctrls:
            letter = control[0]
            outstr += f'\n## {control.capitalize()} Derivatives\n'
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
            outstr += table._repr_markdown_()
        return outstr
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
        if len(self.ctrls) > 0:
            for control in self.ctrls:
                ctrl = self.ctrls[control]
                control = control.capitalize()
                table.add_column(f'{control} (deg)', cfrm, data=[ctrl])
            outstr += table._repr_markdown_()
        if self.gamma is not None:
            table = MDTable()
            table.add_column('Cx', cfrm, data=[self.nfres.Cfrc.x])
            table.add_column('Cy', cfrm, data=[self.nfres.Cfrc.y])
            table.add_column('Cz', cfrm, data=[self.nfres.Cfrc.z])
            outstr += table._repr_markdown_()
            table = MDTable()
            table.add_column('CDi', dfrm, data=[self.nfres.CDi])
            table.add_column('CY', cfrm, data=[self.nfres.CY])
            table.add_column('CL', cfrm, data=[self.nfres.CL])
            table.add_column('e', efrm, data=[self.nfres.e])
            table.add_column('Cl', cfrm, data=[self.nfres.Cl])
            table.add_column('Cm', cfrm, data=[self.nfres.Cm])
            table.add_column('Cn', cfrm, data=[self.nfres.Cn])
            outstr += table._repr_markdown_()
        if self.phi is not None:
            table = MDTable()
            table.add_column('CDi_ff', dfrm, data=[self.trres.CDi])
            table.add_column('CY_ff', cfrm, data=[self.trres.CY])
            table.add_column('CL_ff', cfrm, data=[self.trres.CL])
            table.add_column('e', efrm, data=[self.trres.e])
            table.add_column('Cl_ff', cfrm, data=[self.trres.Cl])
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
    _Cfrc = None
    _Cmom = None
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
            tmp = self.res.sys.amg*self.gamma+self.res.amv
            self._nfmom = elementwise_multiply(self.rhogamma, tmp)
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
    def Cfrc(self):
        if self._Cfrc is None:
            self._Cfrc = self.nffrctot/self.res.qfs/self.res.sys.sref
        return self._Cfrc
    @property
    def Cmom(self):
        if self._Cmom is None:
            self._Cmom = self.nfmomtot/self.res.qfs/self.res.sys.sref
        return self._Cmom
    @property
    def CDi(self):
        if self._CDi is None:
            self._CDi = self.res.acs.dirx*self.Cfrc
        return self._CDi
    @property
    def CY(self):
        if self._CY is None:
            self._CY = self.res.acs.diry*self.Cfrc
        return self._CY
    @property
    def CL(self):
        if self._CL is None:
            self._CL = self.res.acs.dirz*self.Cfrc
        return self._CL
    @property
    def e(self):
        if self._e is None:
            if self.CDi == 0.0:
                self._e = float('nan')
            else:
                from math import pi
                self._e = (self.CL**2+self.CY**2)/pi/self.res.sys.ar/self.CDi
        return self._e
    @property
    def Cl(self):
        if self._Cl is None:
            Cmom = self.res.wcs.vector_to_local(self.Cmom)
            self._Cl = Cmom.x/self.res.sys.bref
        return self._Cl
    @property
    def Cm(self):
        if self._Cm is None:
            Cmom = self.res.wcs.vector_to_local(self.Cmom)
            self._Cm = Cmom.y/self.res.sys.cref
        return self._Cm
    @property
    def Cn(self):
        if self._Cn is None:
            Cmom = self.res.wcs.vector_to_local(self.Cmom)
            self._Cn = Cmom.z/self.res.sys.bref
        return self._Cn

class PhiResult(object):
    res = None
    phi = None
    _trwsh = None
    _trdrg = None
    _trfrc = None
    _trlft = None
    _trmom = None
    _trdrgtot = None
    _trfrctot = None
    _trlfttot = None
    _trmomtot = None
    _CDi = None
    _CY = None
    _CL = None
    _e = None
    _Cl = None
    def __init__(self, res: LatticeResult, phi: matrix):
        self.res = res
        self.phi = phi
    @property
    def trwsh(self):
        if self._trwsh is None:
            self._trwsh = self.res.sys.bvg*self.phi
        return self._trwsh
    @property
    def trdrg(self):
        if self._trdrg is None:
            self._trdrg = self.res.rho*multiply(self.phi, self.res.sys.bdg*self.phi)
        return self._trdrg
    @property
    def trfrc(self):
        if self._trfrc is None:
            self._trfrc = self.res.rho*self.res.speed*multiply(self.phi, self.res.sys.byg)
        return self._trfrc
    @property
    def trlft(self):
        if self._trlft is None:
            self._trlft = self.res.rho*self.res.speed*multiply(self.phi, self.res.sys.blg)
        return self._trlft
    @property
    def trmom(self):
        if self._trmom is None:
            self._trmom = -self.res.rho*self.res.speed*multiply(self.phi, self.res.sys.bmg)
        return self._trmom
    @property
    def trdrgtot(self):
        if self._trdrgtot is None:
            self._trdrgtot = sum(self.trdrg.transpose().tolist()[0])
        return self._trdrgtot
    @property
    def trfrctot(self):
        if self._trfrctot is None:
            self._trfrctot = sum(self.trfrc.transpose().tolist()[0])
        return self._trfrctot
    @property
    def trlfttot(self):
        if self._trlfttot is None:
            self._trlfttot = sum(self.trlft.transpose().tolist()[0])
        return self._trlfttot
    @property
    def trmomtot(self):
        if self._trmomtot is None:
            self._trmomtot = sum(self.trmom.transpose().tolist()[0])
        return self._trmomtot
    @property
    def CDi(self):
        if self._CDi is None:
            self._CDi = self.trdrgtot/self.res.qfs/self.res.sys.sref
        return self._CDi
    @property
    def CY(self):
        if self._CY is None:
            self._CY = self.trfrctot/self.res.qfs/self.res.sys.sref
        return self._CY
    @property
    def CL(self):
        if self._CL is None:
            self._CL = self.trlfttot/self.res.qfs/self.res.sys.sref
        return self._CL
    @property
    def e(self):
        if self._e is None:
            if self.CDi == 0.0:
                self._e = float('nan')
            else:
                from math import pi
                self._e = (self.CL**2+self.CY**2)/pi/self.res.sys.ar/self.CDi
        return self._e
    @property
    def Cl(self):
        if self._Cl is None:
            self._Cl = self.trmomtot/self.res.qfs/self.res.sys.sref/self.res.sys.bref
        return self._Cl
