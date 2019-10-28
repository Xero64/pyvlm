from pyvlm.classes import LatticeSystem
from pyvlm.classes import LatticeTrim

class LevelTrim(object):
    name = None
    sys = None
    gravacc = None
    speed = None
    density = None
    mass = None
    loadfac = None
    _weight = None
    _lift = None
    _dynpres = None
    _acc = None
    _rad = None
    _CL = None
    _prate = None
    _qco2V = None
    def __init__(self, name: str, sys: LatticeSystem):
        self.name = name
        self.sys = sys
        self.gravacc = 9.80665
    def reset(self):
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None
    def set_gravitational_acceleration(self, gravacc: float):
        self.gravacc = gravacc
        self.reset()
    def set_speed_and_density(self, speed: float, density: float):
        self.speed = speed
        self.density = density
        self.reset()
    def set_mass_and_load_factor(self, mass: float, loadfac: float):
        self.mass = mass
        self.loadfac = loadfac
        self.reset()
    def create_trim_result(self):
        ltrm = LatticeTrim(self.name, self.sys)
        ltrm.set_density(rho=self.density)
        ltrm.set_state(speed=self.speed, qco2V=self.qco2V)
        ltrm.set_targets(CLt=self.CL)
        return ltrm
    @property
    def weight(self):
        if self._weight is None:
            self._weight = self.mass*self.gravacc
        return self._weight
    @property
    def lift(self):
        if self._lift is None:
            self._lift = self.loadfac*self.weight
        return self._lift
    @property
    def dynpres(self):
        if self._dynpres is None:
            self._dynpres = self.density*self.speed**2/2
        return self._dynpres
    @property
    def CL(self):
        if self._CL is None:
            self._CL = self.lift/self.dynpres/self.sys.sref
        return self._CL
    @property
    def acc(self):
        if self._acc is None:
            self._acc = (self.loadfac-1)*self.gravacc
        return self._acc
    @property
    def rad(self):
        if self._rad is None:
            if self.acc == 0.0:
                self._rad = float('inf')
            else:
                self._rad = self.speed**2/self.acc
        return self._rad
    @property
    def prate(self):
        if self._prate is None:
            self._prate = self.acc/self.speed
        return self._prate
    @property
    def qco2V(self):
        if self._qco2V is None:
            self._qco2V = self.prate*self.sys.cref/2/self.speed
        return self._qco2V
    def __str__(self):
        from py2md.classes import MDTable
        outstr = '# Trim State '+self.name+' for '+self.sys.name+'\n'
        table = MDTable()
        table.add_column('Speed', '.3f', data=[self.speed])
        table.add_column('Density', '.3f', data=[self.density])
        table.add_column('Dyn. Press.', '.3f', data=[self.dynpres])
        outstr += table._repr_markdown_()
        table = MDTable()
        table.add_column('Mass', '.3f', data=[self.mass])
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
    def _repr_markdown_(self):
        return self.__str__()
