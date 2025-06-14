from time import perf_counter
from typing import TYPE_CHECKING, Any

from numpy import degrees, radians, zeros
from numpy.linalg import norm, solve

from ..tools.mass import Mass
from .latticeresult import LatticeResult as Result

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .latticesystem import LatticeSystem as System

ANGTOL = 30.0


class LatticeTrim(Result):
    targets: dict[str, tuple[str, float]] = None

    def __init__(self, name: str, sys: 'System') -> None:
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
        if display:
            print(f'Cdff = \n{Cdff}\n')
            print(f'Hcur = \n{Hcur}\n')
        Ddff = solve(Hcur, Cdff)
        if display:
            print(f'Ddff = \n{Ddff}\n')
            print(f'Cdff - Hcur@Ddff = \n{Cdff - Hcur@Ddff}\n')
        return Ddff

    def trim(self, crit: float = 1e-6, imax: int = 100, display: bool = False) -> None:
        Ctgt = self.target_Cmat()
        Ccur = self.current_Cmat()
        Cdff = Ctgt - Ccur
        nrmC = norm(Cdff)
        if display:
            print(f'normC = {nrmC}')
        # if display:
        #     outstr = f'{"iter":>4s}  {"d(alpha)":>11s}  {"d(beta)":>11s}'
        #     outstr += f'  {"d(pbo2V)":>11s}  {"d(qco2V)":>11s}  {"d(rbo2V)":>11s}'
        #     for control in self.ctrls:
        #         control = f'd({control})'
        #         outstr += f'  {control:>11s}'
        #     print(outstr)
        count = 0
        while nrmC > crit:
            if display:
                print(f'Iteration {count:d}')
                start = perf_counter()
            Ddff = self.trim_iteration(display = display)
            Dcur = self.current_Dmat() + Ddff
            # if display:
            #     outstr = f'{count+1:4d}'
            #     for i, variable in enumerate(self.targets):
            #         outstr += f'  {Ddff[i]:11.3e}'
            #     print(outstr)
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

    @classmethod
    def from_dict(cls, system: 'System', resdict: dict[str, Any],
                  trim: bool = True) -> 'LatticeTrim':

        from ..tools.trim import (GRAVACC, LevelTrim, LoadTrim, LoopingTrim,
                                TurningTrim)

        name = resdict['name']

        if resdict['trim'] == 'Load Trim':
            trim_condition = LoadTrim(name, system)
            lift = resdict.get('L', None)
            side = resdict.get('Y', None)
            roll = resdict.get('l', None)
            pitch = resdict.get('m', None)
            yaw = resdict.get('n', None)
            trim_condition.set_loads(lift, side, roll, pitch, yaw)

        elif resdict['trim'] == 'Looping Trim':
            trim_condition = LoopingTrim(name, system)
            load_factor = resdict.get('load factor', 1.0)
            trim_condition.set_load_factor(load_factor)

        elif resdict['trim'] == 'Turning Trim':
            trim_condition = TurningTrim(name, system)
            bang = resdict.get('bank angle', 0.0)
            trim_condition.set_bankang(bang)

        elif resdict['trim'] == 'Level Trim':
            trim_condition = LevelTrim(name, system)

        rho = resdict.get('density', 1.0)
        trim_condition.set_density(rho)

        speed = resdict.get('speed', 1.0)
        trim_condition.set_speed(speed)

        gravacc = resdict.get('gravacc', GRAVACC)
        trim_condition.set_gravitational_acceleration(gravacc)

        initstate = {}
        initstate['alpha'] = resdict.get('alpha', 0.0)
        initstate['beta'] = resdict.get('beta', 0.0)
        initstate['pbo2V'] = resdict.get('pbo2V', 0.0)
        initstate['qco2V'] = resdict.get('qco2V', 0.0)
        initstate['rbo2V'] = resdict.get('rbo2V', 0.0)
        trim_condition.set_initial_state(initstate)

        initctrls = {}
        for control in system.ctrls:
            initctrls[control] = resdict.get(control, 0.0)
        trim_condition.set_initial_controls(initctrls)

        mass = resdict.get('mass', None)
        if isinstance(mass, dict):
            mass = Mass(**mass)
        elif isinstance(mass, float):
            mass = Mass(name = trim_condition.name, mass = mass, xcm = system.rref.x,
                        ycm = system.rref.y, zcm = system.rref.z)
        elif mass is None:
            if system.mass is not None:
                mass = system.mass
            else:
                mass = Mass(trim_condition.name, mass = 1.0, xcm = system.rref.x,
                            ycm = system.rref.y, zcm = system.rref.z)
        trim_condition.set_mass(mass)

        trim_result = trim_condition.create_trim_result()

        mach = resdict.get('mach', 0.0)
        trim_result.set_state(mach = mach)

        system.results[name] = trim_result

        if trim:
            trim_result.trim()

        return trim_result
