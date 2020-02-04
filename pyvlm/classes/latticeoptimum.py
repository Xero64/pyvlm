from numpy.matlib import zeros, matrix
from numpy.linalg import solve, norm
from .latticesystem import LatticeSystem
from .latticeresult import LatticeResult
from pygeom.geom3d import Point
from math import degrees
from matplotlib.pyplot import figure

class LatticeOptimum(LatticeResult):
    constr = None
    record = None
    _rhov = None
    _bdg = None
    _adg = None
    def __init__(self, name: str, sys: LatticeSystem):
        super(LatticeOptimum, self).__init__(name, sys)
    def add_constraint(self, param: str, value: float, strplst: list=None, point: Point=None):
        if self.constr is None:
            self.constr = []
        constr = Constraint(self, param, value, point, strplst)
        self.constr.append(constr)
    def add_record(self, param: str, strplst: list=None, point: Point=None):
        if self.record is None:
            self.record = []
        record = Record(self, param, point, strplst)
        self.record.append(record)
    @property
    def rhov(self):
        if self._rhov is None:
            self._rhov = self.rho*self.speed
        return self._rhov
    @property
    def bdg(self):
        return self.sys.bdg
    @property
    def adg(self):
        if self._adg is None:
            self._adg = self.bdg+self.bdg.transpose()
        return self._adg
    def old_iteration(self, pmat: matrix):
        nump = self.sys.nums
        num = nump+len(self.constr)
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
    def optimum_lift_distribution(self, crit=1e-12):
        nump = self.sys.nums
        phi_old = zeros((nump, 1))
        phi, lam = self.old_iteration(phi_old)
        while norm(phi-phi_old) > crit:
            phi_old = phi
            phi, lam = self.old_iteration(phi_old)
        self._phi = phi
        self._lam = lam
        return phi, lam
    def optimum_strip_twist_iteration(self):
        
        nump = len(self.sys.pnls)
        nums = len(self.sys.strps)

        popt = zeros((nums, 1))
        for i, phi in enumerate(self.phi):
            popt[i, 0] = phi

        dafsda = zeros((nump, nums))
        daicda = zeros((nump, nums))

        for i in range(nump):
            for j, strp in enumerate(self.sys.strps):
                for pnl in strp.pnls:
                    k = pnl.lpid
                    if k == i:
                        dafsda[i, j] += self.res.vfs*pnl.dnda
                        daicda[i, j] += self.sys.avc[i, k]*pnl.dnda

        dgda = -solve(self.sys.aic, dafsda+daicda*self.res.phi)

        dpda = zeros((nums, nums))
        for i, strp in enumerate(self.sys.strps):
            for j in range(nums):
                for pnl in strp.pnls:
                    k = pnl.lpid
                    dpda[i, j] += dgda[k, j]

        dal = solve(dpda, popt-self.res.phi)

        for i in range(nums):
            dal[i, 0] = degrees(dal[i, 0])

        da = dal.transpose().tolist()[0]

        al = [strp.angle for i, strp in enumerate(self.sys.strps)]
        alc = [al[i]+da[i] for i in range(nums)]

        self.sys.set_strip_alpha(alc)

        return dal
    def optimum_strip_twist(self, crit=1e-1):
        self.res = LatticeResult(self.name, self.sys)
        self.res.set_state(speed=self.speed)
        self.res.set_density(rho=self.rho)
        dal = self.optimum_strip_twist_iteration()
        nrmdal = norm(dal)
        iteration = 0
        self.res = LatticeResult(self.name, self.sys)
        self.res.set_state(speed=self.speed)
        self.res.set_density(rho=self.rho)
        dal = self.optimum_strip_twist_iteration()
        nrmdal = norm(dal)
        while nrmdal > crit:
            iteration += 1
            print(f'Iteration {iteration} - Convergence {nrmdal}')
            self.res = LatticeResult(self.name, self.sys)
            self.res.set_state(speed=self.speed)
            self.res.set_density(rho=self.rho)
            dal = self.optimum_strip_twist_iteration()
            nrmdal = norm(dal)
        self.res = LatticeResult(self.name, self.sys)
        self.res.set_state(speed=self.speed)
        self.res.set_density(rho=self.rho)
        al = [strp.angle for i, strp in enumerate(self.sys.strps)]
        return al
    def plot_strip_twist_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        y = [strp.pnti.y for i, strp in enumerate(self.sys.strps)]
        al = [strp.angle for i, strp in enumerate(self.sys.strps)]
        ax.plot(y, al, label=f'{self.name} Strip Twist')
        ax.legend()
        return ax
    def return_induced_drag(self):
        return self.rhov*(self.phi.transpose()*self.sys.bdg*self.phi)[0, 0]
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

class Record(object):
    opt = None
    param = None
    point = None
    strplst = None
    _bcv = None
    _bcm = None
    _value = None
    def __init__(self, opt: LatticeOptimum, param: str, pnt: Point=None, strplst: list=None):
        self.opt = opt
        self.param = param
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
                    self._bcv[i, 0] += self.opt.rhov*self.opt.sys.blg[i, 0]
            elif self.param == 'Y':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    self._bcv[i, 0] += self.opt.rhov*self.opt.sys.bdg[i, 0]
            elif self.param == 'l':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    bli = self.opt.sys.blg[i, 0]
                    byi = self.opt.sys.byg[i, 0]
                    ryi = strp.pnti.y-self.point.y
                    rzi = strp.pnti.z-self.point.z
                    self._bcv[i, 0] += self.opt.rhov*(ryi*bli-rzi*byi)
            elif self.param == 'm':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    bli = self.opt.sys.blg[i, 0]
                    rxi = strp.pnti.x-self.point.x
                    self._bcv[i, 0] -= self.opt.rhov*rxi*bli
            elif self.param == 'n':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    byi = self.opt.sys.byg[i, 0]
                    rxi = strp.pnti.x-self.point.x
                    self._bcv[i, 0] += self.opt.rhov*rxi*byi
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
                        self._bcm[i, j] += self.opt.rhov*rzi*self.opt.sys.bdg[i, j]
            elif self.param == 'n':
                self._bcm = zeros((num, num))
                for i in self.strplst:
                    strpi = self.opt.sys.strps[i]
                    ryi = strpi.pnti.y-self.point.y
                    for strpj in self.opt.sys.strps:
                        j = strpj.lsid
                        self._bcm[i, j] -= self.opt.rhov*ryi*self.opt.sys.bdg[i, j]
        return self._bcm
    def return_matrices(self, phi: matrix):
        ach = self.bcv.transpose()
        acv = self.bcv
        if self.bcm is not None:
            ach = ach + phi.transpose()*self.bcm
            acv = acv + (self.bcm.transpose()+self.bcm)*phi
        return acv, ach
    def evaluate(self):
        phi = self.opt.phi
        _, ach = self.return_matrices(phi)
        value = (ach*phi)[0, 0]
        if self._value is None:
            self._value = value
        return value
    @property
    def value(self):
        if self._value is None:
            self._value = self.evaluate()
        return self._value
    def __repr__(self):
        return '<LatticeOptimum Record of {:s}>'.format(self.param)

class Constraint(Record):
    def __init__(self, opt: LatticeOptimum, param: str, value: float, pnt: Point=None, strplst: list=None):
        super(Constraint, self).__init__(opt, param, pnt=pnt, strplst=strplst)
        self._value = value
    def __repr__(self):
        return '<LatticeOptimum Constraint of {:s} to {:}>'.format(self.param, self.value)
