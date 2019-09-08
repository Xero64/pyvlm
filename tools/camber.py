class NACA4(object):
    code = None
    _mt = None
    _mc = None
    _pc = None
    def __init__(self, code: str):
        self.code = code
    @property
    def mt(self):
        if self._mt is None:
            self._mt = float(self.code[2])/10+float(self.code[3])/100
        return self._mt
    @property
    def mc(self):
        if self._mc is None:
            self._mc = float(self.code[0])/100
        return self._mc
    @property
    def pc(self):
        if self._pc is None:
            self._pc = float(self.code[1])/10
        return self._pc
    def return_camber(self, xc: float):
        if xc < self.pc:
            yc = self.mc/self.pc**2*(2*self.pc*xc-xc**2)
        else:
            yc = self.mc/(1-self.pc)**2*((1-2*self.pc)+2*self.pc*xc-xc**2)
        return yc
    def return_camber_slope(self, xc: float):
        if xc < self.pc:
            dydx = self.mc/self.pc**2*(2*self.pc-2*xc)
        else:
            dydx = self.mc/(1-self.pc)**2*(2*self.pc-2*xc)
        return dydx
    def return_camber_angle(self, xc: float):
        from math import degrees, atan
        dydx = self.return_camber_slope(xc)
        return degrees(atan(dydx))
