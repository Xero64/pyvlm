from math import log, pi, copysign

class FlatPlate(object):
    def __init__(self):
        pass
    def return_camber(self, xc: float):
        return 0.0
    def return_camber_slope(self, xc: float):
        return 0.0
    def return_camber_angle(self, xc: float):
        return 0.0
    def __repr__(self):
        return f'<Flat Plate Airfoil>'

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
    def __repr__(self):
        return f'<NACA {self.code:s}>'

class NACA6Series(object):
    code = None
    a = None
    xx = None
    q = None
    cl = None
    h0 = None
    h = None
    g = None
    cff = None
    oma = None
    def __init__(self, code: str):
        if code[0] != '6':
            print('Not a NACA 6 series code.')
            return None
        a = float(code[1])/10
        xx = float(code[-2:])/100
        q = float(code[-3])
        clstr = code[2:-3].lstrip('(').rstrip(')')
        cl = float(clstr)/10
        self.code = code
        self.a = a
        self.xx = xx
        self.q = q
        self.cl = cl
        g = (a**2*(0.5*log(a) - 0.25) + 0.25)/(a - 1)
        self.g = g
        h0 = (0.25 - 0.5*log(1 - a))*(a - 1)
        self.h = h0 + g
        self.cff = cl/(2*pi*(a+1))
    def return_camber(self, xc: float):
        if xc == 0.0:
            yc = 0.0
        else:
            cl = self.cl
            a = self.a
            g = self.g
            h = self.h
            yc = cl*((a - 1)*(g - h*xc - xc*log(xc)) - 0.5*(a - xc)**2*log(abs(a - xc)) + 0.25*(a - xc)**2 + 0.5*(xc - 1)**2*log(1 - xc) - 0.25*(xc - 1)**2)/(2*pi*(a - 1)*(a + 1))
        return yc
    def return_camber_slope(self, xc: float):
        cl = self.cl
        a = self.a
        h = self.h
        cs = cl/(2*pi*(a + 1))
        if xc == 0.0:
            dydx = float('inf')
        elif xc == a:
            dydx = cs*(log(1 - a) - log(a) - h - 1)
        elif xc == 1.0:
            dydx = cs*(copysign((a - xc)**2/2, a - xc) + ((xc-a)/2 - (a - 1)*(h + log(xc) + 1) + (a - xc)*log(abs(a - xc)))*abs(a - xc))/((a - 1)*abs(a - xc))
        else:
            dydx = cs*(copysign((a - xc)**2/2, a - xc) + ((xc-a)/2 - (a - 1)*(h + log(xc) + 1) + (a - xc)*log(abs(a - xc)) + (xc - 1)*log(1 - xc))*abs(a - xc))/((a - 1)*abs(a - xc))
        return dydx
    def return_camber_angle(self, xc: float):
        from math import degrees, atan
        dydx = self.return_camber_slope(xc)
        return degrees(atan(dydx))
    def __repr__(self):
        return f'<NACA {self.code:s}>'
