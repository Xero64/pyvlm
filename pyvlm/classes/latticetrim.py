from time import perf_counter
from typing import TYPE_CHECKING, Any

from numpy import degrees, radians, zeros
from numpy.linalg import norm, solve

from ..tools.mass import Mass
from .latticeresult import LatticeResult

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from .latticesystem import LatticeSystem

ANGTOL = 30.0


class LatticeTrim(LatticeResult):
    targets: dict[str, tuple[str, float]] = None

    def __init__(self, name: str, sys: 'LatticeSystem') -> None:
        super().__init__(name, sys)
        self.targets = {
            'alpha': ('alpha', 0.0),
            'beta': ('beta', 0.0),
            'pbo2V': ('pbo2V', 0.0),
            'qco2V': ('qco2V', 0.0),
            'rbo2V': ('rbo2V', 0.0),
        }
        for control in self.ctrls:
            self.targets[control] = (control, 0.0)

    def set_state(self, mach: float | None = None, speed: float | None = None,
                  alpha: float | None = None, beta: float | None = None,
                  pbo2V: float | None = None, qco2V: float | None = None,
                  rbo2V: float | None = None) -> None:
        super().set_state(mach=mach, speed=speed, alpha=alpha, beta=beta,
                          pbo2V=pbo2V, qco2V=qco2V, rbo2V=rbo2V)
        if self.targets['alpha'][0] == 'alpha':
            self.targets['alpha'] = ('alpha', self.alpha)
        if self.targets['beta'][0] == 'beta':
            self.targets['beta'] = ('beta', self.beta)
        if self.targets['pbo2V'][0] == 'pbo2V':
            self.targets['pbo2V'] = ('pbo2V', self.pbo2V)
        if self.targets['qco2V'][0] == 'qco2V':
            self.targets['qco2V'] = ('qco2V', self.qco2V)
        if self.targets['rbo2V'][0] == 'rbo2V':
            self.targets['rbo2V'] = ('rbo2V', self.rbo2V)

    def set_controls(self, **kwargs: dict[str, float]) -> None:
        super().set_controls(**kwargs)
        for control in kwargs:
            if control in self.targets:
                self.targets[control] = (control, kwargs[control])

    def set_targets(self, CLt: float | None = None, CYt: float | None = None,
                    Clt: float | None = None, Cmt: float | None = None,
                    Cnt: float | None = None) -> None:
        if CLt is not None:
            self.targets['alpha'] = ('CL', CLt)
        if CYt is not None:
            self.targets['beta'] = ('CY', CYt)
        controls = self.ctrls.keys()
        moment = {}
        if Clt is not None:
            moment['Cl'] = Clt
        if Cmt is not None:
            moment['Cm'] = Cmt
        if Cnt is not None:
            moment['Cn'] = Cnt
        for control, (name, value) in zip(controls, moment.items()):
            self.targets[control] = (name, value)

    def set_initial_state(self, initstate: dict[str, float]) -> None:
        if 'alpha' in self.targets:
            if self.targets['alpha'][0] == 'alpha':
                initstate['alpha'] = self.targets['alpha'][1]
        if 'beta' in self.targets:
            if self.targets['beta'][0] == 'beta':
                initstate['beta'] = self.targets['beta'][1]
        if 'pbo2V' in self.targets:
            if self.targets['pbo2V'][0] == 'pbo2V':
                initstate['pbo2V'] = self.targets['pbo2V'][1]
        if 'qco2V' in self.targets:
            if self.targets['qco2V'][0] == 'qco2V':
                initstate['qco2V'] = self.targets['qco2V'][1]
        if 'rbo2V' in self.targets:
            if self.targets['rbo2V'][0] == 'rbo2V':
                initstate['rbo2V'] = self.targets['rbo2V'][1]
        self.initstate = initstate
        super().set_state(**self.initstate)

    def set_initial_controls(self, initctrls: dict[str, float]) -> None:
        for control in self.ctrls:
            if control in self.targets:
                if self.targets[control][0] == control:
                    initctrls[control] = self.targets[control][1]
        self.initctrls = initctrls
        super().set_controls(**self.initctrls)

    def target_Cmat(self) -> 'NDArray':
        Ctgt = zeros(len(self.targets))
        for i, value in enumerate(self.targets.values()):
            Ctgt[i] = value[1]
        return Ctgt

    def current_Cmat(self) -> 'NDArray':
        Ccur = zeros(len(self.targets))
        for i, value in enumerate(self.targets.values()):
            if hasattr(self, value[0]):
                Ccur[i] = getattr(self, value[0])
            elif hasattr(self.nfres, value[0]):
                Ccur[i] = getattr(self.nfres, value[0])
        return Ccur

    def delta_C(self) -> 'NDArray':
        Ctgt = self.target_Cmat()
        Ccur = self.current_Cmat()
        return Ctgt - Ccur

    def current_Dmat(self) -> 'NDArray':
        Dcur = zeros(len(self.targets))
        for i, variable in enumerate(self.targets):
            if variable == 'alpha':
                Dcur[i] = radians(self.alpha)
            elif variable == 'beta':
                Dcur[i] = radians(self.beta)
            elif variable == 'pbo2V':
                Dcur[i] = self.pbo2V
            elif variable == 'qco2V':
                Dcur[i] = self.qco2V
            elif variable == 'rbo2V':
                Dcur[i] = self.rbo2V
            elif variable in self.ctrls:
                Dcur[i] = radians(self.ctrls[variable])
        return Dcur

    def current_Hmat(self) -> 'NDArray':
        num = len(self.targets)
        Hcur = zeros((num, num))
        for i, target in enumerate(self.targets.values()):
            for j, variable in enumerate(self.targets):
                if variable == target[0]:
                    Hcur[i, :] = 0.0
                    Hcur[i, j] = 1.0
                else:
                    if hasattr(self.stres, variable):
                        stvar = getattr(self.stres, variable)
                        if hasattr(stvar, target[0]):
                            Hcur[i, j] = getattr(stvar, target[0])
                    elif variable in self.ctrls:
                        if self.ctrls[variable] < 0.0:
                            ctvar = self.ctresn[variable]
                        else:
                            ctvar = self.ctresp[variable]
                        if hasattr(ctvar, target[0]):
                            Hcur[i, j] = getattr(ctvar, target[0])
        return Hcur

    def trim_iteration(self, display: bool = False) -> 'NDArray':
        Cdff = self.delta_C()
        Hcur = self.current_Hmat()
        Ddff = solve(Hcur, Cdff)
        if display:
            print(f'Cdiff = \n{Cdff}\n')
            print(f'Hcur = \n{Hcur}\n')
            print(f'Ddiff = \n{Ddff}\n')
            print(f'Hcur@Ddiff = \n{Hcur@Ddff}\n')
        return Ddff

    def trim(self, crit: float = 1e-6, imax: int = 100, display: bool = False) -> None:
        Ctgt = self.target_Cmat()
        Ccur = self.current_Cmat()
        Cdff = Ctgt - Ccur
        nrmC = norm(Cdff)
        if display:
            print(f'normC = {nrmC}')
        count = 0
        while nrmC > crit:
            if display:
                print(f'Iteration {count:d}')
                start = perf_counter()
            Ddff = self.trim_iteration(display = display)
            Dcur = self.current_Dmat() + Ddff
            if display:
                finish = perf_counter()
                elapsed = finish - start
                print(f'Trim Internal Iteration Duration = {elapsed:.3f} seconds.')
            state = {}
            ctrls = {}
            for i, variable in enumerate(self.targets):
                if variable == 'alpha':
                    state['alpha'] = degrees(Dcur[i])
                elif variable == 'beta':
                    state['beta'] = degrees(Dcur[i])
                elif variable == 'pbo2V':
                    state['pbo2V'] = Dcur[i]
                elif variable == 'qco2V':
                    state['qco2V'] = Dcur[i]
                elif variable == 'rbo2V':
                    state['rbo2V'] = Dcur[i]
                elif variable in self.ctrls:
                    ctrls[variable] = degrees(Dcur[i])

            failure = False
            for statekey, stateval in state.items():
                if statekey == 'alpha' and abs(stateval) > ANGTOL:
                    print(f'alpha = {stateval:.6f} deg')
                    failure = True
                elif statekey == 'beta' and abs(stateval) > ANGTOL:
                    print(f'beta = {stateval:.6f} deg')
                    failure = True
            for controlkey, controlval in ctrls.items():
                if abs(controlval) > ANGTOL:
                    print(f'{controlkey} = {controlval:.6f} deg')
                    failure = True

            if failure:
                print(f'Convergence failed for {self.name}.')
                return False

            super().set_state(**state)
            super().set_controls(**ctrls)
            Cdff = self.delta_C()
            nrmC = norm(Cdff)
            if display:
                print(f'alpha = {state['alpha']:.6f} deg')
                print(f'beta = {state['beta']:.6f} deg')
                print(f'pbo2V = {state['pbo2V']:.6f}')
                print(f'qco2V = {state['qco2V']:.6f}')
                print(f'rbo2V = {state['rbo2V']:.6f}')
                for control in ctrls:
                    print(f'{control} = {ctrls[control]:.6f} deg')
                print(f'normC = {nrmC}')
            count += 1
            if count >= imax:
                print(f'Convergence failed for {self.name}.')
                return False


def latticetrim_from_dict(system: 'LatticeSystem', resdata: dict[str, Any]) -> LatticeTrim:

    from ..tools.trim import GRAVACC, LevelTrim, LoadTrim, LoopingTrim, TurningTrim

    name = resdata['name']

    if resdata['trim'] == 'Load Trim':
        trim = LoadTrim(name, system)
        lift = resdata.get('L', None)
        side = resdata.get('Y', None)
        roll = resdata.get('l', None)
        pitch = resdata.get('m', None)
        yaw = resdata.get('n', None)
        trim.set_loads(lift, side, roll, pitch, yaw)

    elif resdata['trim'] == 'Looping Trim':
        trim = LoopingTrim(name, system)
        load_factor = resdata.get('load factor', 1.0)
        trim.set_load_factor(load_factor)

    elif resdata['trim'] == 'Turning Trim':
        trim = TurningTrim(name, system)
        bang = resdata.get('bank angle', 0.0)
        trim.set_bankang(bang)

    elif resdata['trim'] == 'Level Trim':
        trim = LevelTrim(name, system)

    rho = resdata.get('density', 1.0)
    trim.set_density(rho)

    speed = resdata.get('speed', 1.0)
    trim.set_speed(speed)

    gravacc = resdata.get('gravacc', GRAVACC)
    trim.set_gravitational_acceleration(gravacc)

    initstate = {}
    initstate['alpha'] = resdata.get('alpha', 0.0)
    initstate['beta'] = resdata.get('beta', 0.0)
    initstate['pbo2V'] = resdata.get('pbo2V', 0.0)
    initstate['qco2V'] = resdata.get('qco2V', 0.0)
    initstate['rbo2V'] = resdata.get('rbo2V', 0.0)
    trim.set_initial_state(initstate)

    initctrls = {}
    for control in system.ctrls:
        initctrls[control] = resdata.get(control, 0.0)
    trim.set_initial_controls(initctrls)

    if isinstance(trim, (LoopingTrim, TurningTrim)):

        mval = 1.0
        xcm, ycm, zcm = system.rref.x, system.rref.y, system.rref.z
        if 'mass' in resdata:
            if isinstance(resdata['mass'], str):
                mass = system.masses[resdata['mass']]
            elif isinstance(resdata['mass'], float):
                if 'rcg' in resdata:
                    rcgdata = resdata['rcg']
                    xcm, ycm, zcm = rcgdata['x'], rcgdata['y'], rcgdata['z']
                mval = resdata['mass']
                mass = Mass(name + ' Mass', mval, xcm, ycm, zcm)
        else:
            mass = Mass(name + ' Mass', mval, xcm, ycm, zcm)

        trim.set_mass(mass)

    trimres = trim.create_trim_result()

    mach = resdata.get('mach', 0.0)
    trimres.set_state(mach = mach)

    system.results[name] = trimres

    trimres.trim()

    return trimres
