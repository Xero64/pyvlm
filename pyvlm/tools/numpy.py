from numpy import divide, pi, reciprocal, zeros
from pygeom.geom3d import Vector

FOURPI = 4.0*pi
TWOPI = 2.0*pi


def numpy_cwdv(pnts: Vector, veca: Vector, vecb: Vector,
               tol: float = 1e-12, betm: float = 1.0) -> Vector:

    rela = pnts - veca
    relb = pnts - vecb

    rela.x = rela.x/betm
    rela.y = rela.y
    rela.z = rela.z
    relb.x = relb.x/betm
    relb.y = relb.y
    relb.z = relb.z

    runa, maga = rela.to_unit(True)
    runb, magb = relb.to_unit(True)

    reca = zeros(maga.shape)
    reciprocal(maga, out=reca, where=maga > tol)

    recb = zeros(magb.shape)
    reciprocal(magb, out=recb, where=magb > tol)

    axb = runa.cross(runb)
    bxc = Vector(0.0, -runb.z, runb.y)
    cxa = Vector(0.0, runa.z, -runa.y)

    abd = 1.0 + runa.dot(runb)
    bcd = 1.0 - runb.x
    cad = 1.0 - runa.x

    facab = zeros(maga.shape)
    divide(reca + recb, abd, where=abd > tol, out=facab)
    velab = axb*facab

    facbc = zeros(maga.shape)
    divide(recb, bcd, where=bcd > tol, out=facbc)
    velbc = bxc*facbc

    facca = zeros(maga.shape)
    divide(reca, cad, where=cad > tol, out=facca)
    velca = cxa*facca

    vel = (velab + velbc + velca)/FOURPI

    return vel
