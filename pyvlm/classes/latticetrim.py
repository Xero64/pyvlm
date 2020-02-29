from time import time
from .latticeresult import LatticeResult, GammaResult
from .latticesystem import LatticeSystem
from numpy.matlib import zeros
from numpy.linalg import norm, inv
from math import degrees, radians
from pygeom.geom3d import Point

class LatticeTrim(LatticeResult):
    CLt = None
    CYt = None
    Clt = None
    Cmt = None
    Cnt = None
    trmfrc = None
    trmmom = None
    trmlft = None
    def __init__(self, name: str, sys: LatticeSystem):
        super(LatticeTrim, self).__init__(name, sys)
        self.set_trim_loads()
    def set_targets(self, CLt: float=0.0, CYt: float=0.0,
                    Clt: float=0.0, Cmt: float=0.0, Cnt: float=0.0):
        self.CLt = CLt
        self.CYt = CYt
        self.Clt = Clt
        self.Cmt = Cmt
        self.Cnt = Cnt
    def set_trim_loads(self, trmfrc: bool=True, trmmom: bool=True, trmlft: bool=False):
        self.trmfrc = trmfrc
        self.trmmom = trmmom
        self.trmlft = trmlft
    def delta_C(self):
        Ctgt = self.target_Cmat()
        Ccur = self.current_Cmat()
        return Ctgt-Ccur
    def target_Cmat(self):
        if self.trmlft:
            Ctgt = zeros((1, 1), dtype=float)
            Ctgt[0, 0] = self.CLt
        elif self.trmfrc and self.trmmom:
            Ctgt = zeros((5, 1), dtype=float)
            Ctgt[0, 0] = self.CLt
            Ctgt[1, 0] = self.CYt
            Ctgt[2, 0] = self.Clt
            Ctgt[3, 0] = self.Cmt
            Ctgt[4, 0] = self.Cnt
        elif self.trmfrc:
            Ctgt = zeros((2, 1), dtype=float)
            Ctgt[0, 0] = self.CLt
            Ctgt[1, 0] = self.CYt
        elif self.trmmom:
            Ctgt = zeros((3, 1), dtype=float)
            Ctgt[0, 0] = self.Clt
            Ctgt[1, 0] = self.Cmt
            Ctgt[2, 0] = self.Cnt
        else:
            Ctgt = zeros((0, 1), dtype=float)
        return Ctgt
    def current_Cmat(self):
        if self.trmlft:
            Ccur = zeros((1, 1), dtype=float)
            Ccur[0, 0] = self.nfres.CL
        elif self.trmfrc and self.trmmom:
            Ccur = zeros((5, 1), dtype=float)
            Ccur[0, 0] = self.nfres.CL
            Ccur[1, 0] = self.nfres.CY
            Ccur[2, 0] = self.nfres.Cl
            Ccur[3, 0] = self.nfres.Cm
            Ccur[4, 0] = self.nfres.Cn
            if self.sys.cdo != 0.0:
                Ccur[0, 0] += self.pdres.CL
                Ccur[1, 0] += self.pdres.CY
                Ccur[2, 0] += self.pdres.Cl
                Ccur[3, 0] += self.pdres.Cm
                Ccur[4, 0] += self.pdres.Cn
        elif self.trmfrc:
            Ccur = zeros((2, 1), dtype=float)
            Ccur[0, 0] = self.nfres.CL
            Ccur[1, 0] = self.nfres.CY
            if self.sys.cdo != 0.0:
                Ccur[0, 0] += self.pdres.CL
                Ccur[1, 0] += self.pdres.CY
        elif self.trmmom:
            Ccur = zeros((3, 1), dtype=float)
            Ccur[0, 0] = self.nfres.Cl
            Ccur[1, 0] = self.nfres.Cm
            Ccur[2, 0] = self.nfres.Cn
            if self.sys.cdo != 0.0:
                Ccur[0, 0] += self.pdres.Cl
                Ccur[1, 0] += self.pdres.Cm
                Ccur[2, 0] += self.pdres.Cn
        else:
            Ccur = zeros((0, 1), dtype=float)
        return Ccur
    def current_Dmat(self):
        if self.trmmom:
            numc = len(self.sys.ctrls)
        else:
            numc = 0
        if self.trmlft:
            Dcur = zeros((1, 1), dtype=float)
            Dcur[0, 0] = radians(self.alpha)
        else:
            num = numc+2
            Dcur = zeros((num, 1), dtype=float)
            Dcur[0, 0] = radians(self.alpha)
            Dcur[1, 0] = radians(self.beta)
        if self.trmmom:
            c = 0
            for control in self.ctrls:
                Dcur[2+c, 0] = radians(self.ctrls[control])
                c += 1
        return Dcur
    def Hmat(self):
        if self.trmmom:
            numc = len(self.sys.ctrls)
        else:
            numc = 0
        num = numc+2
        if self.trmlft:
            H = zeros((1, 1), dtype=float)
            H[0, 0] = self.stres.alpha.CL
        elif self.trmfrc and self.trmmom:
            H = zeros((5, num), dtype=float)
            H[0, 0] = self.stres.alpha.CL
            H[1, 0] = self.stres.alpha.CY
            H[2, 0] = self.stres.alpha.Cl
            H[3, 0] = self.stres.alpha.Cm
            H[4, 0] = self.stres.alpha.Cn
            H[0, 1] = self.stres.beta.CL
            H[1, 1] = self.stres.beta.CY
            H[2, 1] = self.stres.beta.Cl
            H[3, 1] = self.stres.beta.Cm
            H[4, 1] = self.stres.beta.Cn
            ctgam = {}
            c = 0
            for control in self.ctrls:
                if self.ctrls[control] < 0.0:
                    ctgam = self.gctrln_single(control)
                else:
                    ctgam = self.gctrlp_single(control)
                ctres = GammaResult(self, ctgam)
                H[0, 2+c] = ctres.CL
                H[1, 2+c] = ctres.CY
                H[2, 2+c] = ctres.Cl
                H[3, 2+c] = ctres.Cm
                H[4, 2+c] = ctres.Cn
                c += 1
        elif self.trmfrc:
            H = zeros((2, num), dtype=float)
            H[0, 0] = self.stres.alpha.CL
            H[1, 0] = self.stres.alpha.CY
            H[0, 1] = self.stres.beta.CL
            H[1, 1] = self.stres.beta.CY
        elif self.trmmom:
            H = zeros((3, num), dtype=float)
            H[0, 0] = self.stres.alpha.Cl
            H[1, 0] = self.stres.alpha.Cm
            H[2, 0] = self.stres.alpha.Cn
            H[0, 1] = self.stres.beta.Cl
            H[1, 1] = self.stres.beta.Cm
            H[2, 1] = self.stres.beta.Cn
            ctgam = {}
            c = 0
            for control in self.ctrls:
                if self.ctrls[control] < 0.0:
                    ctgam = self.gctrln_single(control)
                else:
                    ctgam = self.gctrlp_single(control)
                ctres = GammaResult(self, ctgam)
                H[0, 2+c] = ctres.Cl
                H[1, 2+c] = ctres.Cm
                H[2, 2+c] = ctres.Cn
                c += 1
        else:
            H = zeros((0, 0), dtype=float)
        return H
    def trim_iteration(self, crit: float=1e-6, imax: int=100):
        Ctgt = self.target_Cmat()
        Ccur = self.current_Cmat()
        Cdff = Ctgt-Ccur
        H = self.Hmat()
        A = H.transpose()*H
        Ainv = inv(A)
        Dcur = self.current_Dmat()
        B = H.transpose()*Cdff
        Ddff = Ainv*B
        Dcur = Dcur+Ddff
        return Dcur
    def trim(self, crit: float=1e-6, imax: int=100, display=False):
        Ctgt = self.target_Cmat()
        Ccur = self.current_Cmat()
        Cdff = Ctgt-Ccur
        nrmC = norm(Cdff)
        if display:
            print(f'normC = {nrmC}')
        i = 0
        while nrmC > crit:
            if display:
                print(f'Iteration {i:d}')
                start = time()
            Dcur = self.trim_iteration(crit, imax)
            if display:
                finish = time()
                elapsed = finish-start
                print(f'Trim Internal Iteration Duration = {elapsed:.3f} seconds.')
            if Dcur is False:
                return
            alpha = degrees(Dcur[0, 0])
            if self.trmlft:
                beta = self.beta
            else:
                beta = degrees(Dcur[1, 0])
            if self.trmmom:
                ctrls = {}
                c = 0
                for control in self.ctrls:
                    ctrls[control] = degrees(Dcur[2+c, 0])
                    c += 1
                self.set_controls(**ctrls)
            self.set_state(alpha=alpha, beta=beta)
            Ccur = self.current_Cmat()
            Cdff = Ctgt-Ccur
            nrmC = norm(Cdff)
            if display:
                print(f'alpha = {alpha:.6f} deg')
                print(f'beta = {beta:.6f} deg')
                for control in ctrls:
                    print(f'{control} = {ctrls[control]:.6f} deg')
                print(f'normC = {nrmC}')
            i += 1
            if i >= imax:
                print(f'Convergence failed for {self.name}.')
                return False

def latticetrim_from_json(lsys: LatticeSystem, resdata: dict):
    from pyvlm.tools.trim import LoopingTrim, TurningTrim
    from pyvlm.tools.mass import Mass
    name = resdata['name']

    if resdata['trim'] == 'Looping Trim':
        trim = LoopingTrim(name, lsys)
        n = 1.0
        if 'load factor' in resdata:
            n = resdata['load factor']
        trim.set_load_factor(n)
    elif resdata['trim'] == 'Turning Trim':
        trim = TurningTrim(name, lsys)
        bang = 0.0
        if 'bank angle' in resdata:
            bang = resdata['bank angle']
        trim.set_bank_angle(bang)

    rho = 1.0
    if 'density' in resdata:
        rho = resdata['density']
    speed = 1.0
    if 'speed' in resdata:
        speed = resdata['speed']
    trim.set_speed_and_density(speed, rho)

    if 'gravacc' in resdata:
        g = resdata['gravacc']
        trim.set_gravitational_acceleration(g)

    m = 1.0
    xcm, ycm, zcm = lsys.rref.x, lsys.rref.y, lsys.rref.z
    if 'mass' in resdata:
        if isinstance(resdata['mass'], str):
            mass = lsys.masses[resdata['mass']]
        elif isinstance(resdata['mass'], float):
            if 'rcg' in resdata:
                rcgdata = resdata['rcg']
                xcm, ycm, zcm = rcgdata['x'], rcgdata['y'], rcgdata['z']
            m = resdata['mass']
            mass = Mass(name + ' Mass', m, xcm, ycm, zcm)
    else:
        mass = Mass(name + ' Mass', m, xcm, ycm, zcm)
    
    trim.set_mass(mass)

    ltrm = trim.create_trim_result()

    lsys.results[name] = ltrm

    ltrm.trim()

    return ltrm
