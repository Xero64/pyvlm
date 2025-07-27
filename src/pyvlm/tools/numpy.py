from numpy import divide, pi, reciprocal, zeros
from pygeom.geom3d import Vector

FOURPI = 4.0*pi


def numpy_cwdv(pnts: Vector, veca: Vector, vecb: Vector, *,
               tol: float = 1e-12, betm: float = 1.0) -> Vector:

    a = pnts - veca
    b = pnts - vecb

    a.x = a.x/betm
    a.y = a.y
    a.z = a.z
    b.x = b.x/betm
    b.y = b.y
    b.z = b.z

    ah, am = a.to_unit(True)
    bh, bm = b.to_unit(True)

    ar = zeros(am.shape)
    reciprocal(am, out=ar, where=am > tol)

    br = zeros(bm.shape)
    reciprocal(bm, out=br, where=bm > tol)

    axb = ah.cross(bh)
    bxc = Vector(0.0, -bh.z, bh.y)
    cxa = Vector(0.0, ah.z, -ah.y)

    abd = 1.0 + ah.dot(bh)
    bcd = 1.0 - bh.x
    cad = 1.0 - ah.x

    facab = zeros(am.shape)
    divide(ar + br, abd, where=abd > tol, out=facab)

    facbc = zeros(am.shape)
    divide(br, bcd, where=bcd > tol, out=facbc)

    facca = zeros(am.shape)
    divide(ar, cad, where=cad > tol, out=facca)

    vel = (axb*facab + bxc*facbc + cxa*facca)/FOURPI

    return vel
