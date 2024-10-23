from typing import TYPE_CHECKING

from pygeom.geom3d import Vector

from ..classes import LatticeTrim
from .mass import Mass, MassCollection

if TYPE_CHECKING:
    from ..classes import LatticeSystem


GRAVACC = 9.80665

class LoopingTrim():
    name: str = None
    sys: 'LatticeSystem' = None
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

    def __init__(self, name: str, sys: 'LatticeSystem') -> None:
        self.name = name
        self.sys = sys
        self.gravacc = GRAVACC

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    def set_gravacc(self, value: float) -> None:
        self.gravacc = value
        self.reset()

    def set_speed(self, value: float) -> None:
        self.speed = value
        self.reset()

    def set_density(self, value: float) -> None:
        self.density = value
        self.reset()

    def set_loadfac(self, value: float) -> None:
        self.loadfac = value
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

    def create_trim_result(self) -> LatticeTrim:
        ltrm = LatticeTrim(self.name, self.sys)
        ltrm.set_density(rho=self.density)
        ltrm.set_state(speed=self.speed, qco2V=self.qco2V)
        ltrm.set_targets(CLt=self.CL)
        rcg = Vector(self.mass.xcm, self.mass.ycm, self.mass.zcm)
        ltrm.set_cg(rcg)
        return ltrm

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

    def __str__(self) -> str:
        from py2md.classes import MDTable
        outstr = '# Looping Trim State '+self.name+' for '+self.sys.name+'\n'
        table = MDTable()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Mass', '.3f', data=[self.mass.mass])
        table.add_column('Grav. Acc.', '.5f', data=[self.gravacc])
        table.add_column('Weight', '.3f', data=[self.weight])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Load Factor', '.3f', data=[self.loadfac])
        table.add_column('Lift', '.3f', data=[self.lift])
        table.add_column('CL', '.5f', data=[self.CL])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Acceleration', '.3f', data=[self.acc])
        table.add_column('Radius', '.3f', data=[self.rad])
        table.add_column('Pitch Rate', '.5f', data=[self.prate])
        outstr += table._repr_markdown_()
        return outstr

    def _repr_markdown_(self) -> str:
        return self.__str__()


class TurningTrim():
    name: str = None
    sys: 'LatticeSystem' = None
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

    def __init__(self, name: str, sys: 'LatticeSystem') -> None:
        self.name = name
        self.sys = sys
        self.gravacc = GRAVACC

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    def set_gravacc(self, value: float) -> None:
        self.gravacc = value
        self.reset()

    def set_speed(self, value: float) -> None:
        self.speed = value
        self.reset()

    def set_density(self, value: float) -> None:
        self.density = value
        self.reset()

    def set_bankang(self, value: float) -> None:
        self.bankang = value
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

    def create_trim_result(self) -> LatticeTrim:
        ltrm = LatticeTrim(self.name, self.sys)
        ltrm.set_density(rho=self.density)
        ltrm.set_state(speed=self.speed, qco2V=self.qco2V, rbo2V=self.rbo2V)
        ltrm.set_targets(CLt=self.CL)
        rcg = Vector(self.mass.xcm, self.mass.ycm, self.mass.zcm)
        ltrm.set_cg(rcg)
        return ltrm

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
                fac = (self.loadfac**2-1.0)/self.loadfac
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

    def __str__(self) -> str:
        from py2md.classes import MDTable
        outstr = '# Turning Trim State '+self.name+' for '+self.sys.name+'\n'
        table = MDTable()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Mass', '.3f', data=[self.mass.mass])
        table.add_column('Grav. Acc.', '.5f', data=[self.gravacc])
        table.add_column('Weight', '.3f', data=[self.weight])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Bank Angle (deg)', '.1f', data=[self.bankang])
        table.add_column('Load Factor', '.3f', data=[self.loadfac])
        table.add_column('Lift', '.3f', data=[self.lift])
        table.add_column('CL', '.5f', data=[self.CL])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Acceleration', '.3f', data=[self.acc])
        table.add_column('Turn Radius', '.3f', data=[self.rad])
        table.add_column('Pitch Rate', '.5f', data=[self.prate])
        table.add_column('Roll Rate', '.5f', data=[self.rrate])
        outstr += table._repr_markdown_()
        return outstr

    def _repr_markdown_(self) -> str:
        return self.__str__()


class LevelTrim():
    name: str = None
    sys: 'LatticeSystem' = None
    gravacc: float = None
    speed: float = None
    density: float = None
    mass: 'Mass | MassCollection' = None
    _weight: float = None
    _lift: float = None
    _dynpres: float = None
    _CL: float = None

    def __init__(self, name: str, sys: 'LatticeSystem') -> None:
        self.name = name
        self.sys = sys
        self.gravacc = GRAVACC

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    def set_gravacc(self, value: float) -> None:
        self.gravacc = value
        self.reset()

    def set_speed(self, value: float) -> None:
        self.speed = value
        self.reset()

    def set_density(self, value: float) -> None:
        self.density = value
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

    def create_trim_result(self) -> 'LatticeTrim':
        lres = LatticeTrim(self.name, self.sys)
        lres.set_density(rho=self.density)
        lres.set_state(speed=self.speed)
        return lres

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

    def trim_speed_from_CL(self, CL: float) -> float:
        if self.mass is not None and self.density is not None:
            W = self.weight
            S = self.sys.sref
            rho = self.density
            self.speed = (W/S/rho/CL*2)**0.5
            self._CL = CL

    def __str__(self) -> str:
        from py2md.classes import MDTable
        outstr = '# Level Trim State '+self.name+' for '+self.sys.name+'\n'
        table = MDTable()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Mass', '.3f', data=[self.mass.mass])
        table.add_column('Grav. Acc.', '.5f', data=[self.gravacc])
        table.add_column('Weight', '.3f', data=[self.weight])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Lift', '.3f', data=[self.lift])
        table.add_column('CL', '.5f', data=[self.CL])
        outstr += table._repr_markdown_()
        return outstr

    def _repr_markdown_(self) -> str:
        return self.__str__()
