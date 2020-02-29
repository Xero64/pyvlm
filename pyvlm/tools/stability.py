from cmath import sqrt

g = 9.80665 # m/s*82

class StabilityApproximation(object):
    res = None
    mass = None
    def __init__(self, res, mass):
        self.res = res
        self.mass = mass
    def phugoid(self):
        Xu = self.res.stres.u.nffrctot.x/self.mass.mass
        Zu = self.res.stres.u.nffrctot.z/self.mass.mass
        u0 = self.res.vfs.x
        return (Xu+sqrt(Xu**2+4*g*Zu/u0))/2, (Xu-sqrt(Xu**2+4*g*Zu/u0))/2
    def short_period(self):
        lmnq = self.res.wcs.vector_to_local(self.res.stres.q.nfmomtot)
        Mq = lmnq.y/self.mass.Iyy
        lmnw = self.res.wcs.vector_to_local(self.res.stres.w.nfmomtot)
        Mw = lmnw.y/self.mass.Iyy
        u0 = self.res.vfs.x
        return (Mq+sqrt(Mq**2+4*u0*Mw))/2, (Mq-sqrt(Mq**2+4*u0*Mw))/2
    def roll_subsidence(self):
        lmnp = self.res.wcs.vector_to_local(self.res.stres.p.nfmomtot)
        Lp = lmnp.x/self.mass.Ixx
        return Lp
    def spiral(self):
        lmnr = self.res.wcs.vector_to_local(self.res.stres.r.nfmomtot)
        Lr = lmnr.x/self.mass.Ixx
        Nr = lmnr.z/self.mass.Izz
        lmnv = self.res.wcs.vector_to_local(self.res.stres.v.nfmomtot)
        Lv = lmnv.x/self.mass.Ixx
        Nv = lmnv.z/self.mass.Izz
        return Nr-Nv*Lr/Lv
    def dutch_roll(self):
        xyzr = self.res.stres.r.nffrctot
        Yr = xyzr.y/self.mass.mass
        lmnr = self.res.wcs.vector_to_local(self.res.stres.r.nfmomtot)
        Nr = lmnr.z/self.mass.Izz
        xyzv = self.res.stres.v.nffrctot
        Yv = xyzv.y/self.mass.mass
        lmnv = self.res.wcs.vector_to_local(self.res.stres.v.nfmomtot)
        Nv = lmnv.z/self.mass.Izz
        u0 = self.res.vfs.x
        return Nr+Yv+sqrt(u0*Nv+Yv*Nr-Nv*Yr), Nr+Yv-sqrt(u0*Nv+Yv*Nr-Nv*Yr)
