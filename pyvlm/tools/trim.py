from typing import TYPE_CHECKING

from py2md.classes import MDReport
from pygeom.geom3d import Vector

from ..classes import LatticeTrim as Trim
from .mass import Mass, MassCollection

if TYPE_CHECKING:
    from ..classes import LatticeSystem as System


GRAVACC = 9.80665

class LoopingTrim():
    name: str = None
    sys: 'System' = None
    gravacc: float = None
    speed: float = None
    density: float = None
    mass: 'Mass | MassCollection' = None
    loadfac: float = None
    _weight: float = None
    _lift: float = None
    _dynpres: float = None
    _acc: float = None
    _rad: float = None
    _CL: float = None
    _prate: float = None
    _qco2V: float = None
    initstate: dict[str, float] = None
    initctrls: dict[str, float] = None

    def __init__(self, name: str, sys: 'System') -> None:
        self.name = name
        self.sys = sys
        self.gravacc = GRAVACC
        self.initstate = {}
        self.initctrls = {}

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    def set_gravitational_acceleration(self, gravacc: float) -> None:
        self.gravacc = gravacc
        self.reset()

    def set_speed(self, speed: float) -> None:
        self.speed = speed
        self.reset()

    def set_density(self, value: float) -> None:
        self.density = value
        self.reset()

    def set_load_factor(self, loadfac: float) -> None:
        self.loadfac = loadfac
        self.reset()

    def set_mass(self, mass: 'Mass | MassCollection | str | float') -> None:
        if isinstance(mass, str):
            self.mass = self.sys.masses[mass]
        elif isinstance(mass, float):
            self.mass = Mass(self.name + ' Mass', mass,
                             self.sys.rref.x, self.sys.rref.y, self.sys.rref.z)
        elif isinstance(mass, (Mass, MassCollection)):
            self.mass = mass
        self.reset()

    def set_initial_state(self, initstate: dict[str, float]) -> None:
        initstate.setdefault('alpha', 0.0)
        initstate.setdefault('beta', 0.0)
        initstate.setdefault('pbo2V', 0.0)
        initstate.setdefault('qco2V', 0.0)
        initstate.setdefault('rbo2V', 0.0)
        self.initstate = initstate

    def set_initial_controls(self, initctrls: dict[str, float]) -> None:
        for control in self.sys.ctrls:
            initctrls.setdefault(control, 0.0)
        self.initctrls = initctrls

    def create_trim_result(self) -> Trim:
        trim = Trim(self.name, self.sys)
        trim.set_density(rho=self.density)
        trim.set_state(speed=self.speed, qco2V=self.qco2V)
        trim.set_targets(CLt = self.CL, CYt = 0.0, Clt = 0.0, Cmt = 0.0, Cnt = 0.0)
        rcg = Vector(self.mass.xcm, self.mass.ycm, self.mass.zcm)
        trim.set_cg(rcg)
        trim.set_initial_state(self.initstate)
        trim.set_initial_controls(self.initctrls)
        return trim

    @property
    def weight(self) -> float:
        if self._weight is None:
            self._weight = self.mass.mass*self.gravacc
        return self._weight

    @property
    def lift(self) -> float:
        if self._lift is None:
            self._lift = self.loadfac*self.weight
        return self._lift

    @property
    def dynpres(self) -> float:
        if self._dynpres is None:
            self._dynpres = self.density*self.speed**2/2
        return self._dynpres

    @property
    def CL(self) -> float:
        if self._CL is None:
            self._CL = self.lift/self.dynpres/self.sys.sref
        return self._CL

    @property
    def acc(self) -> float:
        if self._acc is None:
            self._acc = (self.loadfac - 1.0)*self.gravacc
        return self._acc

    @property
    def rad(self) -> float:
        if self._rad is None:
            if self.acc == 0.0:
                self._rad = float('inf')
            else:
                self._rad = self.speed**2/self.acc
        return self._rad

    @property
    def prate(self) -> float:
        if self._prate is None:
            self._prate = self.acc/self.speed
        return self._prate

    @property
    def qco2V(self) -> float:
        if self._qco2V is None:
            self._qco2V = self.prate*self.sys.cref/2/self.speed
        return self._qco2V

    def to_mdobj(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        table = report.add_table()
        table.add_column('Mass', '.3f', data=[self.mass.mass])
        table.add_column('Grav. Acc.', '.5f', data=[self.gravacc])
        table.add_column('Weight', '.3f', data=[self.weight])
        table = report.add_table()
        table.add_column('Load Factor', '.3f', data=[self.loadfac])
        table.add_column('Lift', '.3f', data=[self.lift])
        table.add_column('CL', '.5f', data=[self.CL])
        table = report.add_table()
        table.add_column('Acceleration', '.3f', data=[self.acc])
        table.add_column('Radius', '.3f', data=[self.rad])
        table.add_column('Pitch Rate', '.5f', data=[self.prate])
        return report

    def __str__(self) -> str:
        return self.to_mdobj().__str__()

    def _repr_markdown_(self) -> str:
        return self.to_mdobj()._repr_markdown_()


class TurningTrim():
    name: str = None
    sys: 'System' = None
    gravacc: float = None
    speed: float = None
    density: float = None
    mass: 'Mass | MassCollection' = None
    bankang: float = None
    _loadfac: float = None
    _weight: float = None
    _lift: float = None
    _dynpres: float = None
    _acc: float = None
    _rad: float = None
    _CL: float = None
    _prate: float = None
    _rrate: float = None
    _qco2V: float = None
    _rbo2V: float = None
    initstate: dict[str, float] = None
    initctrls: dict[str, float] = None

    def __init__(self, name: str, sys: 'System') -> None:
        self.name = name
        self.sys = sys
        self.gravacc = GRAVACC
        self.initstate = {}
        self.initctrls = {}

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    def set_gravitational_acceleration(self, gravacc: float):
        self.gravacc = gravacc
        self.reset()

    def set_speed(self, speed: float) -> None:
        self.speed = speed
        self.reset()

    def set_density(self, density: float) -> None:
        self.density = density
        self.reset()

    def set_bankang(self, bankang: float) -> None:
        self.bankang = bankang
        self.reset()

    def set_mass(self, mass: 'Mass | MassCollection | str | float') -> None:
        if isinstance(mass, str):
            self.mass = self.sys.masses[mass]
        elif isinstance(mass, float):
            self.mass = Mass(self.name + ' Mass', mass,
                             self.sys.rref.x, self.sys.rref.y, self.sys.rref.z)
        elif isinstance(mass, (Mass, MassCollection)):
            self.mass = mass
        self.reset()

    def set_initial_state(self, initstate: dict[str, float]) -> None:
        initstate.setdefault('alpha', 0.0)
        initstate.setdefault('beta', 0.0)
        initstate.setdefault('pbo2V', 0.0)
        initstate.setdefault('qco2V', 0.0)
        initstate.setdefault('rbo2V', 0.0)
        self.initstate = initstate

    def set_initial_controls(self, initctrls: dict[str, float]) -> None:
        for control in self.sys.ctrls:
            initctrls.setdefault(control, 0.0)
        self.initctrls = initctrls

    def create_trim_result(self) -> Trim:
        trim = Trim(self.name, self.sys)
        trim.set_density(rho=self.density)
        trim.set_state(speed = self.speed, qco2V = self.qco2V, rbo2V = self.rbo2V)
        trim.set_targets(CLt = self.CL, CYt = 0.0, Clt = 0.0, Cmt = 0.0, Cnt = 0.0)
        rcg = Vector(self.mass.xcm, self.mass.ycm, self.mass.zcm)
        trim.set_cg(rcg)
        trim.set_initial_state(self.initstate)
        trim.set_initial_controls(self.initctrls)
        return trim

    @property
    def loadfac(self) -> float:
        if self._loadfac is None:
            from math import cos, radians
            brad = radians(self.bankang)
            self._loadfac = 1.0/cos(brad)
        return self._loadfac

    @property
    def weight(self) -> float:
        if self._weight is None:
            self._weight = self.mass.mass*self.gravacc
        return self._weight

    @property
    def lift(self) -> float:
        if self._lift is None:
            self._lift = self.loadfac*self.weight
        return self._lift

    @property
    def dynpres(self) -> float:
        if self._dynpres is None:
            self._dynpres = self.density*self.speed**2/2
        return self._dynpres

    @property
    def CL(self) -> float:
        if self._CL is None:
            self._CL = self.lift/self.dynpres/self.sys.sref
        return self._CL

    @property
    def acc(self) -> float:
        if self._acc is None:
            self._acc = (self.loadfac**2 - 1.0)**0.5*self.gravacc
        return self._acc

    @property
    def rad(self) -> float:
        if self._rad is None:
            if self.acc == 0.0:
                self._rad = float('inf')
            else:
                self._rad = self.speed**2/self.acc
        return self._rad

    @property
    def prate(self) -> float:
        if self._prate is None:
            if self.acc != 0.0:
                fac = (self.loadfac**2 - 1.0)/self.loadfac
                self._prate = self.gravacc/self.speed*fac
            else:
                self._prate = 0.0
        return self._prate

    @property
    def rrate(self) -> float:
        if self._rrate is None:
            if self.acc != 0.0:
                self._rrate = self.acc/self.speed/self.loadfac
            else:
                self._rrate = 0.0
        return self._rrate

    @property
    def qco2V(self) -> float:
        if self._qco2V is None:
            self._qco2V = self.prate*self.sys.cref/2/self.speed
        return self._qco2V

    @property
    def rbo2V(self) -> float:
        if self._rbo2V is None:
            self._rbo2V = self.rrate*self.sys.bref/2/self.speed
        return self._rbo2V

    def to_mdobj(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        table = report.add_table()
        table.add_column('Mass', '.3f', data=[self.mass.mass])
        table.add_column('Grav. Acc.', '.5f', data=[self.gravacc])
        table.add_column('Weight', '.3f', data=[self.weight])
        table = report.add_table()
        table.add_column('Bank Angle (deg)', '.1f', data=[self.bankang])
        table.add_column('Load Factor', '.3f', data=[self.loadfac])
        table.add_column('Lift', '.3f', data=[self.lift])
        table.add_column('CL', '.5f', data=[self.CL])
        table = report.add_table()
        table.add_column('Acceleration', '.3f', data=[self.acc])
        table.add_column('Turn Radius', '.3f', data=[self.rad])
        table.add_column('Pitch Rate', '.5f', data=[self.prate])
        table.add_column('Roll Rate', '.5f', data=[self.rrate])
        return report

    def __str__(self) -> str:
        return self.to_mdobj().__str__()

    def _repr_markdown_(self) -> str:
        return self.to_mdobj()._repr_markdown_()


class LevelTrim():
    name: str = None
    sys: 'System' = None
    gravacc: float = None
    speed: float = None
    density: float = None
    mass: 'Mass | MassCollection' = None
    _weight: float = None
    _lift: float = None
    _dynpres: float = None
    _CL: float = None
    initstate: dict[str, float] = None
    initctrls: dict[str, float] = None

    def __init__(self, name: str, sys: 'System') -> None:
        self.name = name
        self.sys = sys
        self.gravacc = GRAVACC
        self.initstate = {}
        self.initctrls = {}

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    def set_gravitational_acceleration(self, gravacc: float) -> None:
        self.gravacc = gravacc
        self.reset()

    def set_density(self, density: float) -> None:
        self.density = density
        self.reset()

    def set_speed(self, speed: float):
        self.speed = speed
        self.reset()

    def set_mass(self, mass: 'Mass | MassCollection | str | float') -> None:
        if isinstance(mass, str):
            self.mass = self.sys.masses[mass]
        elif isinstance(mass, float):
            self.mass = Mass(self.name + ' Mass', mass,
                             self.sys.rref.x, self.sys.rref.y, self.sys.rref.z)
        elif isinstance(mass, (Mass, MassCollection)):
            self.mass = mass
        self.reset()

    def set_initial_state(self, initstate: dict[str, float]) -> None:
        initstate.setdefault('alpha', 0.0)
        initstate.setdefault('beta', 0.0)
        initstate.setdefault('pbo2V', 0.0)
        initstate.setdefault('qco2V', 0.0)
        initstate.setdefault('rbo2V', 0.0)
        self.initstate = initstate

    def set_initial_controls(self, initctrls: dict[str, float]) -> None:
        for control in self.sys.ctrls:
            initctrls.setdefault(control, 0.0)
        self.initctrls = initctrls

    def create_trim_result(self) -> Trim:
        trim = Trim(self.name, self.sys)
        trim.set_density(rho=self.density)
        trim.set_state(speed=self.speed)
        trim.set_targets(CLt = self.CL)
        rcg = Vector(self.mass.xcm, self.mass.ycm, self.mass.zcm)
        trim.set_cg(rcg)
        trim.set_initial_state(self.initstate)
        trim.set_initial_controls(self.initctrls)
        return trim

    @property
    def weight(self) -> float:
        if self._weight is None:
            self._weight = self.mass.mass*self.gravacc
        return self._weight

    @property
    def lift(self) -> float:
        if self._lift is None:
            self._lift = self.weight
        return self._lift

    @property
    def dynpres(self) -> float:
        if self._dynpres is None:
            self._dynpres = self.density*self.speed**2/2
        return self._dynpres

    @property
    def CL(self) -> float:
        if self._CL is None:
            self._CL = self.lift/self.dynpres/self.sys.sref
        return self._CL

    def trim_speed_from_CL(self, CL: float) -> None:
        if self.mass is not None and self.density is not None:
            W = self.weight
            S = self.sys.sref
            rho = self.density
            self.speed = (W/S/rho/CL*2)**0.5
            self._CL = CL

    def to_mdobj(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        table = report.add_table()
        table.add_column('Mass', '.3f', data=[self.mass.mass])
        table.add_column('Grav. Acc.', '.5f', data=[self.gravacc])
        table.add_column('Weight', '.3f', data=[self.weight])
        table = report.add_table()
        table.add_column('Lift', '.3f', data=[self.lift])
        table.add_column('CL', '.5f', data=[self.CL])
        return report

    def __str__(self) -> str:
        return self.to_mdobj().__str__()

    def _repr_markdown_(self) -> str:
        return self.to_mdobj()._repr_markdown_()


class LoadTrim():
    name: str = None
    sys: 'System' = None
    speed: float = None
    density: float = None
    L: float = None
    Y: float = None
    l: float = None
    m: float = None
    n: float = None
    initstate: dict[str, float] = None
    initctrls: dict[str, float] = None
    _dynpres: float = None
    _CL: float = None
    _CY: float = None
    _Cl: float = None
    _Cm: float = None
    _Cn: float = None

    def __init__(self, name: str, sys: 'System') -> None:
        self.name = name
        self.sys = sys
        self.initstate = {}
        self.initctrls = {}

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    def set_density(self, density: float) -> None:
        self.density = density
        self.reset()

    def set_speed(self, speed: float):
        self.speed = speed
        self.reset()

    def set_loads(self, L: float = None, Y: float = None,
                  l: float = None, m: float = None, n: float = None) -> None:
        self.L = L
        self.Y = Y
        self.l = l
        self.m = m
        self.n = n

    def set_initial_state(self, initstate: dict[str, float]) -> None:
        initstate.setdefault('alpha', 0.0)
        initstate.setdefault('beta', 0.0)
        initstate.setdefault('pbo2V', 0.0)
        initstate.setdefault('qco2V', 0.0)
        initstate.setdefault('rbo2V', 0.0)
        self.initstate = initstate

    def set_initial_controls(self, initctrls: dict[str, float]) -> None:
        for control in self.sys.ctrls:
            initctrls.setdefault(control, 0.0)
        self.initctrls = initctrls

    def create_trim_result(self) -> Trim:
        trim = Trim(self.name, self.sys)
        trim.set_density(rho=self.density)
        trim.set_state(speed=self.speed)
        trim.set_targets(CLt = self.CL, CYt = self.CY,
                         Clt = self.Cl, Cmt = self.Cm, Cnt = self.Cn)
        trim.set_initial_state(self.initstate)
        trim.set_initial_controls(self.initctrls)
        return trim

    @property
    def dynpres(self) -> float:
        if self._dynpres is None:
            self._dynpres = self.density*self.speed**2/2
        return self._dynpres

    @property
    def CL(self) -> float:
        if self._CL is None:
            if self.L is not None:
                self._CL = self.L/self.dynpres/self.sys.sref
        return self._CL

    @property
    def CY(self) -> float:
        if self._CY is None:
            if self.Y is not None:
                self._CY = self.Y/self.dynpres/self.sys.sref
        return self._CY

    @property
    def Cl(self) -> float:
        if self._Cl is None:
            if self.l is not None:
                self._Cl = self.l/self.dynpres/self.sys.sref/self.sys.bref
        return self._Cl

    @property
    def Cm(self) -> float:
        if self._Cm is None:
            if self.m is not None:
                self._Cm = self.m/self.dynpres/self.sys.sref/self.sys.cref
        return self._Cm

    @property
    def Cn(self) -> float:
        if self._Cn is None:
            if self.n is not None:
                self._Cn = self.n/self.dynpres/self.sys.sref/self.sys.bref
        return self._Cn

    def to_mdobj(self) -> MDReport:
        report = MDReport()
        table = report.add_table()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        table = report.add_table()
        table.add_column('L', '.3f', data=[self.L])
        table.add_column('Y', '.3f', data=[self.Y])
        table.add_column('l', '.3f', data=[self.l])
        table.add_column('m', '.3f', data=[self.m])
        table.add_column('n', '.3f', data=[self.n])
        table = report.add_table()
        table.add_column('CL', '.5f', data=[self.CL])
        table.add_column('CY', '.5f', data=[self.CY])
        table.add_column('Cl', '.5f', data=[self.Cl])
        table.add_column('Cm', '.5f', data=[self.Cm])
        table.add_column('Cn', '.5f', data=[self.Cn])
        return report

    def __str__(self) -> str:
        return self.to_mdobj().__str__()

    def _repr_markdown_(self) -> str:
        return self.to_mdobj()._repr_markdown_()
