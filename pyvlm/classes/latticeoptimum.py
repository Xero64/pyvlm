from typing import List

from matplotlib.pyplot import figure
from numpy import asarray, degrees
from numpy.linalg import norm, solve
from numpy.matlib import matrix, zeros
from pygeom.geom3d import Vector

from .latticeresult import LatticeResult


class LatticeOptimum(LatticeResult):
    constr = None
    record = None
    phiopt = None
    lamopt = None
    _bdg = None
    _adg = None

    def add_constraint(self, param: str, value: float, strplst: list=None, point: Vector=None):
        if self.constr is None:
            self.constr = []
        constr = Constraint(self, param, value, point, strplst)
        self.constr.append(constr)

    def add_record(self, param: str, strplst: list=None, point: Vector=None):
        if self.record is None:
            self.record = []
        record = Record(self, param, point, strplst)
        self.record.append(record)

    def set_target_phi(self, phi: List[float]) -> None:
        if len(phi) != len(self.sys.strps):
            raise Exception('The length of phi must equal the number of strips.')
        self.phiopt = asarray(phi)
        self._phi = self.phiopt

    def set_target_lift_force_distribution(self, ltgt: list, rho: float, speed: float, mach: float=0.0):
        if len(ltgt) != len(self.sys.strps):
            raise Exception('The length of l must equal the number of strips.')
        phitgt = [li/rho/speed for li in ltgt]
        self.set_density(rho)
        self.set_state(mach=mach, speed=speed)
        self.set_target_phi(phitgt)

    @property
    def bdg(self):
        return self.sys.bdg

    @property
    def adg(self):
        if self._adg is None:
            self._adg = self.bdg + self.bdg.transpose()
        return self._adg

    def old_iteration(self, pmat: matrix):
        nump = self.sys.nums
        num = nump + len(self.constr)
        amat = zeros((num, num))
        bmat = zeros((num, 1))
        amat[0:nump, 0:nump] = self.adg
        for i, constr in enumerate(self.constr):
            acv, ach = constr.return_matrices(pmat)
            amat[0:nump, nump+i] = acv
            amat[nump+i, 0:nump] = ach
            bmat[nump+i, 0] = constr.value
        xmat = solve(amat, bmat)
        phi = xmat[0:nump, 0]
        lam = xmat[nump:num, 0]
        return phi, lam

    def optimum_lift_force_distribution(self, crit=1e-12):
        nump = self.sys.nums
        phi_old = zeros((nump, 1))
        phi, lam = self.old_iteration(phi_old)
        while norm(phi-phi_old) > crit:
            phi_old = phi
            phi, lam = self.old_iteration(phi_old)
        self._phi = phi
        self.phiopt = phi
        self.lamopt = lam
        return phi, lam

    def optimum_strip_twist_iteration(self):

        nump = len(self.sys.pnls)
        nums = len(self.sys.strps)

        dafsda = zeros((nump, nums))
        daicda = zeros((nump, nums))

        aic = self.sys.aic(self.mach)
        avc = self.sys.avc(self.mach)

        for i in range(nump):
            for j, strp in enumerate(self.sys.strps):
                for pnl in strp.pnls:
                    k = pnl.lpid
                    if k == i:
                        dafsda[i, j] += self.vfs*pnl.dnda
                        daicda[i, j] += avc[i, k]*pnl.dnda

        dgda = -solve(aic, dafsda + daicda*self.phi)

        dpda = zeros((nums, nums))
        for i, strp in enumerate(self.sys.strps):
            for j in range(nums):
                for pnl in strp.pnls:
                    k = pnl.lpid
                    dpda[i, j] += dgda[k, j]

        dphi = self.phiopt-self.phi

        dal = solve(dpda, dphi)

        for i in range(nums):
            dal[i, 0] = degrees(dal[i, 0])

        da = dal.transpose().tolist()[0]

        al = [strp.twist for i, strp in enumerate(self.sys.strps)]
        alc = [al[i]+da[i] for i in range(nums)]

        self.sys.set_strip_alpha(alc)
        self._ungam = None
        self._gamma = None
        self._phi = None
        self._trres = None

        return dal

    def optimum_strip_twist(self, crit=1e-1):
        self._ungam = None
        self._gamma = None
        self._phi = None
        self._trres = None
        dal = self.optimum_strip_twist_iteration()
        nrmdal = norm(dal)
        iteration = 0
        dal = self.optimum_strip_twist_iteration()
        nrmdal = norm(dal)
        while nrmdal > crit:
            iteration += 1
            print(f'Iteration {iteration} - Convergence alg {nrmdal}')
            dal = self.optimum_strip_twist_iteration()
            # self.reset()
            nrmdal = norm(dal)
        al = [strp.twist for i, strp in enumerate(self.sys.strps)]
        return al

    def plot_target_phi_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.sys.srfcs[0].strpy, self.phi, label=self.name+' Target')
        ax.legend()
        return ax

    def plot_strip_twist_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        y = [strp.pnti.y for i, strp in enumerate(self.sys.strps)]
        al = [strp.twist for i, strp in enumerate(self.sys.strps)]
        ax.plot(y, al, label=f'{self.name} Strip Twist')
        ax.legend()
        return ax

    def return_induced_drag(self):
        temp = self.phi.transpose()*self.sys.bdg*self.phi
        return self.rho*self.speed*temp[0, 0]

    def __repr__(self):
        return '<LatticeOptimum {:s}>'.format(self.name)

    def __str__(self):
        from py2md.classes import MDTable
        outstr = '# '+self.name+'\n'
        table = MDTable()
        table.add_column('Speed', 'g', data=[self.speed])
        table.add_column('Density', 'g', data=[self.rho])
        outstr += table._repr_markdown_()
        if self._phi is not None:
            table = MDTable()
            table.add_column('Label', 's')
            table.add_column('Type', 's')
            table.add_column('Value', 'g')
            table.add_row(['Di', 'Objective', self.return_induced_drag()])
            if self.constr is not None:
                for constr in self.constr:
                    val = constr.evaluate()
                    table.add_row([constr.param, 'Constraint', val])
            if self.record is not None:
                for record in self.record:
                    val = record.evaluate()
                    table.add_row([record.param, 'Record', val])
            outstr += table._repr_markdown_()
        return outstr

    def _repr_markdown_(self):
        return self.__str__()

class Record():
    opt = None
    param = None
    point = None
    strplst = None
    _bcv = None
    _bcm = None
    _rhov = None
    _value = None

    def __init__(self, opt: LatticeOptimum, param: str, pnt: Vector=None, strplst: list=None):
        self.opt = opt
        self.param = param
        self.point = pnt
        self.strplst = strplst
        self.update()

    def update(self):
        if self.strplst is None:
            self.strplst = [strp.lsid for strp in self.opt.sys.strps]
        if self.point is None:
            self.point = self.opt.sys.rref

    @property
    def bcv(self):
        if self._bcv is None:
            num = self.opt.sys.nums
            if self.param == 'L':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    self._bcv[i, 0] += self.rhov*self.opt.sys.blg[i, 0]
            elif self.param == 'Y':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    self._bcv[i, 0] += self.rhov*self.opt.sys.bdg[i, 0]
            elif self.param == 'l':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    bli = self.opt.sys.blg[i, 0]
                    byi = self.opt.sys.byg[i, 0]
                    ryi = strp.pnti.y-self.point.y
                    rzi = strp.pnti.z-self.point.z
                    self._bcv[i, 0] += self.rhov*(ryi*bli-rzi*byi)
            elif self.param == 'm':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    bli = self.opt.sys.blg[i, 0]
                    rxi = strp.pnti.x-self.point.x
                    self._bcv[i, 0] -= self.rhov*rxi*bli
            elif self.param == 'n':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    byi = self.opt.sys.byg[i, 0]
                    rxi = strp.pnti.x-self.point.x
                    self._bcv[i, 0] += self.rhov*rxi*byi
        return self._bcv

    @property
    def bcm(self):
        if self._bcm is None:
            num = self.opt.sys.nums
            if self.param == 'm':
                self._bcm = zeros((num, num))
                for i in self.strplst:
                    strpi = self.opt.sys.strps[i]
                    rzi = strpi.pnti.z-self.point.z
                    for strpj in self.opt.sys.strps:
                        j = strpj.lsid
                        self._bcm[i, j] += self.rhov*rzi*self.opt.sys.bdg[i, j]
            elif self.param == 'n':
                self._bcm = zeros((num, num))
                for i in self.strplst:
                    strpi = self.opt.sys.strps[i]
                    ryi = strpi.pnti.y-self.point.y
                    for strpj in self.opt.sys.strps:
                        j = strpj.lsid
                        self._bcm[i, j] -= self.rhov*ryi*self.opt.sys.bdg[i, j]
        return self._bcm

    @property
    def rhov(self):
        if self._rhov is None:
            self._rhov = self.opt.rho*self.opt.speed
        return self._rhov

    def return_matrices(self, phi: matrix):
        ach = self.bcv.transpose()
        acv = self.bcv
        if self.bcm is not None:
            ach = ach + phi.transpose()*self.bcm
            acv = acv + (self.bcm.transpose() + self.bcm)*phi
        return acv, ach

    def evaluate(self) -> float:
        phi = self.opt.phi
        _, ach = self.return_matrices(phi)
        value = (ach*phi)[0, 0]
        if self._value is None:
            self._value = value
        return value

    @property
    def value(self) -> float:
        if self._value is None:
            self._value = self.evaluate()
        return self._value

    def __repr__(self):
        return '<LatticeOptimum Record of {:s}>'.format(self.param)

class Constraint(Record):

    def __init__(self, opt: LatticeOptimum, param: str, value: float,
                 pnt: Vector=None, strplst: list=None) -> None:
        super(Constraint, self).__init__(opt, param, pnt=pnt, strplst=strplst)
        self._value = value

    def __repr__(self) -> str:
        return '<LatticeOptimum Constraint of {:s} to {:}>'.format(self.param, self.value)
