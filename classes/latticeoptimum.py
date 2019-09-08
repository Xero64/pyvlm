from numpy.matlib import zeros, matrix
from numpy.linalg import solve, norm
from .latticesystem import LatticeSystem
from .latticeresult import LatticeResult
from pygeom.geom3d import Point, Vector, zero_vector
from math import degrees
from matplotlib.pyplot import figure

class LatticeOptimum(object):
    name = None
    sys = None
    sym = None
    ind = None
    constr = None
    record = None
    rho = None
    speed = None
    res = None
    _phi = None
    _pmat = None
    _rhov = None
    _bdg = None
    _adg = None
    def __init__(self, name: str, sys: LatticeSystem, sym: bool=True):
        self.name = name
        self.sys = sys
        self.sym = sym
        self.update()
    def update(self):
        self.ind = []
        for strp in self.sys.strps:
            self.ind.append(None)
        if self.sym:
            ind = 0
            for strp in self.sys.strps:
                if strp.msid is not None:
                    self.ind[strp.lsid] = ind
                    self.ind[strp.msid] = ind
                    ind += 1
            for strp in self.sys.strps:
                if self.ind[strp.lsid] is None:
                    if strp.msid is None:
                        self.ind[strp.lsid] = ind
                        ind += 1
        else:
            ind = 0
            for strp in self.sys.strps:
                self.ind[strp.lsid] = ind
                ind += 1
    def set_conditions(self, speed: float=1.0, rho: float=1.0):
        self.speed = speed
        self.rho = rho
    def add_constraint(self, param: str, value: float, strplst: list=None, point: Point=None):
        if self.constr is None:
            self.constr = []
        constr = Constraint(self, param, value, point, strplst)
        self.constr.append(constr)
    def add_record(self, param: str, strplst: list=None, point: Point=None):
        if self.record is None:
            self.record = []
        record = Constraint(self, param, None, point, strplst)
        self.record.append(record)
    def initial_pmat(self):
        nump = max(self.ind)+1
        return zeros((nump, 1))
    @property
    def rhov(self):
        if self._rhov is None:
            self._rhov = self.rho*self.speed
        return self._rhov
    @property
    def bdg(self):
        if self._bdg is None:
            num = max(self.ind)+1
            self._bdg = zeros((num, num))
            for i, indi in enumerate(self.ind):
                for j, indj in enumerate(self.ind):
                    self._bdg[indi, indj] += self.rhov*self.sys.bdg[i, j]
        return self._bdg
    @property
    def adg(self):
        if self._adg is None:
            self._adg = self.bdg+self.bdg.transpose()
        return self._adg
    def old_iteration(self, pmat: matrix):
        nump = max(self.ind)+1
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
        pmat = xmat[0:nump, 0]
        lmat = xmat[nump:num, 0]
        return pmat, lmat
    def optimum_lift_distribution(self, crit=1e-12):
        pmat_cur = self.initial_pmat()
        pmat_new, lmat = self.old_iteration(pmat_cur)
        while norm(pmat_new-pmat_cur) > crit:
            pmat_cur = pmat_new
            pmat_new, lmat = self.old_iteration(pmat_cur)
        phi = [pmat_new[ind, 0] for ind in self.ind]
        lam = [lmat[i, 0] for i in range(len(lmat))]
        self._pmat = pmat_new
        self._lmat = lmat
        return phi, lam
    @property
    def phi(self):
        if self._phi is None:
            self._phi, _ = self.optimum_lift_distribution()
        return self._phi
    @property
    def pmat(self):
        if self._pmat is None:
            nump = max(self.ind)+1
            self._pmat = zeros((nump, 1))
            for i, ind in enumerate(self.ind):
                self._pmat[ind, 0] = self.phi[i]
        return self._pmat
    def set_phi(self, phi: list):
        if len(phi) != len(self.sys.strps):
            raise Exception('The length of phi must equal the number of strips.')
        self._phi = phi
    def set_lift_distribution(self, l: list, rho: float, speed: float):
        if len(l) != len(self.sys.strps):
            raise Exception('The length of l must equal the number of strips.')
        self._phi = [li/rho/speed for li in l]
        self.speed = speed
        self.rho = rho
    def return_induced_drag(self):
        return (self.pmat.transpose()*self.bdg*self.pmat)[0, 0]
    def print_report(self):
        print(f'{self.name} Report')
        Di = self.return_induced_drag()
        print(f'Di = {Di:.6f}')
        if self.constr is not None:
            for constr in self.constr:
                val = constr.evaluate(self.pmat)
                print(f'{constr.param} = {val:.6f}')
        if self.record is not None:
            for record in self.record:
                val = record.evaluate(self.pmat)
                print(f'{record.param} = {val:.6f}')
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

        dgda = -solve(self.sys.aic, dafsda+daicda*self.res.pmat)

        dpda = zeros((nums, nums))
        for i, strp in enumerate(self.sys.strps):
            for j in range(nums):
                for pnl in strp.pnls:
                    k = pnl.lpid
                    dpda[i, j] += dgda[k, j]

        dal = solve(dpda, popt-self.res.pmat)

        for i in range(nums):
            dal[i, 0] = degrees(dal[i, 0])

        da = dal.transpose().tolist()[0]

        al = [strp.angle for i, strp in enumerate(self.sys.strps)]
        alc = [al[i]+da[i] for i in range(nums)]

        self.sys.set_strip_alpha(alc)

        return dal
    def optimum_strip_twist(self, crit=1e-1):
        self.res = LatticeResult(self.name, self.sys)
        self.res.set_conditions(speed=self.speed, rho=self.rho)
        dal = self.optimum_strip_twist_iteration()
        nrmdal = norm(dal)
        iteration = 0
        self.res = LatticeResult(self.name, self.sys)
        self.res.set_conditions(speed=self.speed, rho=self.rho)
        dal = self.optimum_strip_twist_iteration()
        nrmdal = norm(dal)
        while nrmdal > crit:
            iteration += 1
            print(f'Iteration {iteration} - Convergence {nrmdal}')
            self.res = LatticeResult(self.name, self.sys)
            self.res.set_conditions(speed=self.speed, rho=self.rho)
            dal = self.optimum_strip_twist_iteration()
            nrmdal = norm(dal)
        self.res = LatticeResult(self.name, self.sys)
        self.res.set_conditions(speed=self.speed, rho=self.rho)
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
    def __repr__(self):
        return '<LatticeOptimum {:s}>'.format(self.name)
    def __str__(self):
        from py2md.classes import MDTable
        outstr = '# '+self.name+'\n'
        table = MDTable()
        table.add_column('Speed [m/s]', 'g', data=[self.speed])
        table.add_column('Rho [kg/m**3]', 'g', data=[self.rho])
        outstr += table._repr_markdown_()
        if self._phi is not None:
            table = MDTable()
            table.add_column('Label', 's')
            table.add_column('Type', 's')
            table.add_column('Value', 'g')
            table.add_row(['Di', 'Objective', self.return_induced_drag()])
            if self.constr is not None:
                for constr in self.constr:
                    val = constr.evaluate(self.pmat)
                    table.add_row([constr.param, 'Constraint', val])
            if self.record is not None:
                for record in self.record:
                    val = record.evaluate(self.pmat)
                    table.add_row([record.param, 'Record', val])
            outstr += table._repr_markdown_()
        return outstr
    def _repr_markdown_(self):
        return self.__str__()

class Constraint(object):
    opt = None
    param = None
    value = None
    point = None
    strplst = None
    _bcv = None
    _bcm = None
    def __init__(self, opt: LatticeOptimum, param: str, value: float, pnt: Point=None, strplst: list=None):
        self.opt = opt
        self.param = param
        self.value = value
        self.strplst = strplst
        self.update()
    def update(self):
        if self.strplst is None:
            self.strplst = []
            for strp in self.opt.sys.strps:
                self.strplst.append(strp.lsid)
        elif self.strplst == 'Mirrored':
            self.strplst = []
            for strp in self.opt.sys.strps:
                if strp.msid is not None:
                    self.strplst.append(strp.lsid)
        if self.point is None:
            self.point = self.opt.sys.rref
    @property
    def bcv(self):
        if self._bcv is None:
            num = max(self.opt.ind)+1
            if self.param == 'L':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    ind = self.opt.ind[i]
                    self._bcv[ind, 0] += self.opt.rhov*self.opt.sys.blg[i, 0]
            elif self.param == 'Y':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    ind = self.opt.ind[i]
                    self._bcv[ind, 0] += self.opt.rhov*self.opt.sys.bdg[i, 0]
            elif self.param == 'l':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    ind = self.opt.ind[i]
                    bli = self.opt.sys.blg[i, 0]
                    byi = self.opt.sys.byg[i, 0]
                    ryi = strp.pnti.y-self.point.y
                    rzi = strp.pnti.z-self.point.z
                    self._bcv[ind, 0] += self.opt.rhov*(ryi*bli-rzi*byi)
            elif self.param == 'm':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    ind = self.opt.ind[i]
                    bli = self.opt.sys.blg[i, 0]
                    rxi = strp.pnti.x-self.point.x
                    self._bcv[ind, 0] -= self.opt.rhov*rxi*bli
            elif self.param == 'n':
                self._bcv = zeros((num, 1))
                for i in self.strplst:
                    strp = self.opt.sys.strps[i]
                    ind = self.opt.ind[i]
                    byi = self.opt.sys.byg[i, 0]
                    rxi = strp.pnti.x-self.point.x
                    self._bcv[ind, 0] += self.opt.rhov*rxi*byi
        return self._bcv
    @property
    def bcm(self):
        if self._bcm is None:
            num = max(self.opt.ind)+1
            if self.param == 'm':
                self._bcm = zeros((num, num))
                for i in self.strplst:
                    strpi = self.opt.sys.strps[i]
                    indi = self.opt.ind[i]
                    rzi = strpi.pnti.z-self.point.z
                    for j, indj in enumerate(self.opt.ind):
                        bdi = self.opt.sys.bdg[i, j]
                        self._bcm[indi, indj] += self.opt.rhov*rzi*bdi
            elif self.param == 'n':
                self._bcm = zeros((num, num))
                for i in self.strplst:
                    strpi = self.opt.sys.strps[i]
                    indi = self.opt.ind[i]
                    ryi = strpi.pnti.y-self.point.y
                    for j, indj in enumerate(self.opt.ind):
                        bdi = self.opt.sys.bdg[i, j]
                        self._bcm[indi, indj] -= self.opt.rhov*ryi*bdi
        return self._bcm
    def return_matrices(self, pmat: matrix):
        ach = self.bcv.transpose()
        acv = self.bcv
        if self.bcm is not None:
            ach = ach + pmat.transpose()*self.bcm
            acv = acv + (self.bcm.transpose()+self.bcm)*pmat
        return acv, ach
    def evaluate(self, pmat: matrix):
        _, ach = self.return_matrices(pmat)
        value = (ach*pmat)[0, 0]
        return value
