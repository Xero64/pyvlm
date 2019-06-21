"""Module containing function for 1D interpolation"""

class Function(object):
    x = None
    y = None
    m = None
    dy2 = None
    def __init__(self, x, y):
        if isinstance(x, list) and isinstance(y, list):
            lenx = len(x)
            leny = len(y)
            if lenx == leny and lenx > 1 and leny > 1:
                ascending = True
                for i in range(lenx-1):
                    if x[i] > x[i+1]:
                        ascending = False
                if ascending:
                    self.x = x
                    self.y = y
        self.calc_d2y()
    def calc_m(self):
        num = len(self.x)
        self.m = []
        for i in range(num-1):
            dx = (self.x[i+1]-self.x[i])
            if dx == 0.:
                self.m.append(0.)
            else:
                dy = (self.y[i+1]-self.y[i])
                self.m.append(dy/dx)
    def linear_interp(self, xi):
        if self.m == None:
            self.calc_m()
        num = len(self.x)
        if isinstance(xi, list):
            numi = len(xi)
            yi = [None]*numi
            for i in range(numi):
                for j in range(num-1):
                    if xi[i] >= self.x[j] and xi[i] <= self.x[j+1]:
                        yi[i] = self.y[j]+self.m[j]*(xi[i]-self.x[j])
                        break
            return yi
        if isinstance(xi, int):
            xi = float(xi)
        if isinstance(xi, float):
            for j in range(num-1):
                if xi >= self.x[j] and xi <= self.x[j+1]:
                    yi = self.y[j]+self.m[j]*(xi-self.x[j])
                    break
            return yi
    def cumtrapz(self):
        A = [0.]
        for i in range(1,len(self.x)):
            A.append(A[i-1]+(self.y[i]+self.y[i-1])/2*(self.x[i]-self.x[i-1]))
        return A
    def calc_d2y(self):
        num = len(self.x)
        a = [0. for i in range(num)]
        b = [1. for i in range(num)]
        c = [0. for i in range(num)]
        r = [0. for i in range(num)]
        for i in range(1, num-1):
            dxA = self.x[i]-self.x[i-1]
            dxB = self.x[i+1]-self.x[i]
            dyA = self.y[i]-self.y[i-1]
            dyB = self.y[i+1]-self.y[i]
            a[i] = dxA/6
            b[i] = (dxA+dxB)/3
            c[i] = dxB/6
            r[i] = dyB/dxB-dyA/dxA
        gm = [0. for i in range(num)]
        bt = b[0]
        self.d2y = [0. for i in range(num)]
        self.d2y[0] = r[0]/bt
        for i in range(1, num):
            gm[i] = c[i-1]/bt
            bt = b[i]-a[i]*gm[i]
            self.d2y[i] = (r[i]-a[i]*self.d2y[i-1])/bt
        for i in range(num-2, 0, -1):
            self.d2y[i] -= gm[i+1]*self.d2y[i+1]
    def cubic_interp(self, xi):
        num = len(self.x)
        if isinstance(xi, list):
            numi = len(xi)
            yi = [None]*numi
            for i in range(numi):
                for j in range(num-1):
                    if xi[i] >= self.x[j] and xi[i] <= self.x[j+1]:
                        dx = (self.x[j+1]-self.x[j])
                        A = (self.x[j+1]-xi[i])/dx
                        B = (xi[i]-self.x[j])/dx
                        C = (A**3-A)*dx**2/6
                        D = (B**3-B)*dx**2/6
                        yi[i] = A*self.y[j]+B*self.y[j+1]+C*self.d2y[j]+D*self.d2y[j-1]
            return yi
        if isinstance(xi, int):
            xi = float(xi)
        if isinstance(xi, float):
            for j in range(num-1):
                if xi >= self.x[j] and xi <= self.x[j+1]:
                    dx = (self.x[j+1]-self.x[j])
                    A = (self.x[j+1]-xi)/dx
                    B = (xi-self.x[j])/dx
                    C = (A**3-A)*dx**2/6
                    D = (B**3-B)*dx**2/6
                    yi = A*self.y[j]+B*self.y[j+1]+C*self.d2y[j]+D*self.d2y[j-1]
            return yi
    def cubic_interp_eval(self, xi=None, level=0, sv1=0., sv2=0.):
        if xi == None:
            if level == 0:
                return self.y
        num = len(self.x)
        if xi == None:
            if level == -2:
                return self.d2y
        iy0 = [0. for i in range(num)]
        i2y0 = [0. for i in range(num)]
        iy0[0] = sv1
        i2y0[0] = sv2
        for i in range(1, num):
            xa = self.x[i-1]
            xb = self.x[i]
            ya = self.y[i-1]
            yb = self.y[i]
            ua = self.d2y[i-1]
            ub = self.d2y[i]
            iy0[i] = iy0[i-1]
            iy0[i] += (xb-xa)/2*ya
            iy0[i] += (xb-xa)/2*yb
            iy0[i] += (xa-xb)**3/24*ua
            iy0[i] += (xa-xb)**3/24*ub
            i2y0[i] = iy0[i-1]*(xb-xa)+i2y0[i-1]
            i2y0[i] += (2*xb**2+2*xa*xb-xa**2)/6*ya
            i2y0[i] += (xb**2-2*xa*xb-2*xa**2)/6*yb
            i2y0[i] += (7*xa**4-28*xa**3*xb+12*xa**2*xb**2+32*xa*xb**3-8*xb**4)/360*ua
            i2y0[i] += (8*xa**4-32*xa**3*xb-12*xa**2*xb**2+28*xa*xb**3-7*xb**4)/360*ub
        if xi == None:
            if level == 1:
                return iy0
            if level == 2:
                return i2y0
        return iy0, i2y0

def function_from_dict(dinp):
    if isinstance(dinp, dict):
        x = [k for k in sorted(dinp)]
        y = [dinp[k] for k in sorted(dinp)]
        return Function(x, y)
