from pygeom.geom3d import Vector


class LatticeGrid(Vector):
    gid: int = None

    def __init__(self, x: float, y: float, z: float) -> None:
        super().__init__(x, y, z)
