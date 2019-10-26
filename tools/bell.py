from math import pi, sqrt, acos, sin
from matplotlib.pyplot import figure

class Bell(object):
    lift = None
    span = None
    y = None
    _ym = None
    _s = None
    _sm = None
    _th = None
    _thm = None
    _speed = None
    _rho = None
    def __init__(self, span: float, y: list):
        self.span = span
        self.y = y
    def set_lift(self, lift: float):
        self.lift = lift
    def set_speed(self, speed: float):
        self._speed = speed
    def set_density(self, rho: float):
        self._rho = rho
    def set_ym(self, ym: list):
        self._ym = ym
    @property
    def s(self):
        if self._s is None:
            self._s = [2*yi/self.span for yi in self.y]
        return self._s
    @property
    def th(self):
        if self._th is None:
            self._th = [acos(si) for si in self.s]
        return self._th
    @property
    def speed(self):
        if self._speed is None:
            self._speed = 1.0
        return self._speed
    @property
    def rho(self):
        if self._rho is None:
            self._rho = 1.0
        return self._rho
    @property
    def drag(self):
        L = self.lift
        b = self.span
        V = self.speed
        rho = self.rho
        return 8*L**2/(3*pi*V**2*b**2*rho)
    @property
    def ym(self):
        if self._ym is None:
            self._ym = []
            for i in range(len(self.y)-1):
                a, b = i, i+1
                ya = self.y[a]
                yb = self.y[b]
                self._ym.append((ya+yb)/2)
        return self._ym
    @property
    def sm(self):
        if self._sm is None:
            self._sm = [2*yi/self.span for yi in self.ym]
        return self._sm
    @property
    def thm(self):
        if self._thm is None:
            self._thm = [acos(si) for si in self.sm]
        return self._thm
    def return_lift_forces(self):
        L = self.lift
        Li = []
        for i in range(len(self.th)-1):
            a, b = i, i+1
            tha = self.th[a]
            thb = self.th[b]
            Li.append(L*(3*tha - 3*thb - 2*sin(2*tha) + sin(4*tha)/4 + 2*sin(2*thb) - sin(4*thb)/4)/(3*pi))
        return Li
    def return_delta_y(self):
        dy = []
        for i in range(len(self.y)-1):
            a, b = i, i+1
            ya = self.y[a]
            yb = self.y[b]
            dy.append(yb-ya)
        return dy
    def lift_distribution(self):
        L = self.lift
        b = self.span
        return [16*L*sin(th)**3/(3*pi*b) for th in self.th]
    def drag_distribution(self):
        L = self.lift
        b = self.span
        V = self.speed
        rho = self.rho
        return [32*L**2*(sin(th) - sin(3*th))*sin(th)**2/(3*pi**2*V**2*b**3*rho) for th in self.th]
    def wash_distribution(self):
        L = self.lift
        b = self.span
        V = self.speed
        rho = self.rho
        return [8*L*(2*s**2-1)/(pi*V*b**2*rho) for s in self.s]
    def trefftz_wash_distribution(self):
        L = self.lift
        b = self.span
        V = self.speed
        rho = self.rho
        return [4*L*(2*s**2-1)/(pi*V*b**2*rho) for s in self.s]
    def shear_force_distribution(self):
        L = self.lift
        th = [thi if thi <= pi/2 else thi-pi for thi in self.th]
        return [L*(3*thi - 2*sin(2*thi) + sin(4*thi)/4)/(3*pi) for thi in th]
    def bending_moment_distribution(self):
        L = self.lift
        b = self.span
        th = [thi if thi <= pi/2 else thi-pi for thi in self.th]
        return [L*b*(3*thi**2 - sin(thi)**4 - 3*sin(thi)**2)*sin(thi)/(12*pi) for thi in th]
    def root_shear_force(self):
        L = self.lift
        return L/2
    def root_bending_moment(self):
        L = self.lift
        b = self.span
        return L*b*(-16 + 3*pi**2)/(48*pi)
    def return_phi(self):
        L = self.lift
        b = self.span
        V = self.speed
        rho = self.rho
        return [16*L*sin(th)**3/(3*pi*b)/rho/V for th in self.thm]
    def __str__(self):
        sf = self.root_shear_force()
        bm = self.root_bending_moment()
        outstr = ''
        outstr += f'Span = {self.span:g}\n'
        outstr += f'Lift = {self.lift:g}\n'
        outstr += f'Drag = {self.drag:g}\n'
        outstr += f'Speed = {self.speed:g}\n'
        outstr += f'Density = {self.rho:g}\n'
        outstr += f'Root Shear Force = {sf:g}\n'
        outstr += f'Root Bending Moment = {bm:g}\n'
        return outstr
    def plot_lift_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.y, self.lift_distribution(), label='Bell Theory')
        ax.legend()
        return ax
    def plot_drag_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.y, self.drag_distribution(), label='Bell Theory')
        ax.legend()
        return ax
    def plot_wash_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.y, self.wash_distribution(), label='Bell Theory')
        ax.legend()
        return ax
    def plot_shear_force_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.y, self.shear_force_distribution(), label='Bell Theory')
        ax.legend()
        return ax
    def plot_bending_moment_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.y, self.bending_moment_distribution(), label='Bell Theory')
        ax.legend()
        return ax
    def plot_trefftz_wash_distribution(self, ax=None):
        if ax is None:
            fig = figure(figsize=(12, 8))
            ax = fig.gca()
            ax.grid(True)
        ax.plot(self.y, self.trefftz_wash_distribution(), label='Bell Theory')
        ax.legend()
        return ax
