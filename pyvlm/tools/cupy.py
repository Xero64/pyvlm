from time import perf_counter

from pygeom.geom3d import Vector

try:
    from cupy import ElementwiseKernel, asarray

    cupy_cwdv_kernel = ElementwiseKernel(
        '''float64 px, float64 py, float64 pz,
        float64 ax, float64 ay, float64 az,
        float64 bx, float64 by, float64 bz,
        float64 tol, float64 betm''',
        '''float64 vx, float64 vy, float64 vz''',
        '''

        const double fourpi = 12.566370614359172463991854;

        double rax = (px - ax) / betm;
        double ray = py - ay;
        double raz = pz - az;
        double ram = sqrt(rax*rax + ray*ray + raz*raz);

        double hax = 0.0;
        double hay = 0.0;
        double haz = 0.0;
        double rar = 0.0;
        if (ram > tol) {
            rar = 1.0 / ram;
            hax = rax / ram;
            hay = ray / ram;
            haz = raz / ram;
        }

        double rbx = (px - bx) / betm;
        double rby = py - by;
        double rbz = pz - bz;
        double rbm = sqrt(rbx*rbx + rby*rby + rbz*rbz);

        double hbx = 0.0;
        double hby = 0.0;
        double hbz = 0.0;
        double rbr = 0.0;
        if (rbm > tol) {
            rbr = 1.0 / rbm;
            hbx = rbx / rbm;
            hby = rby / rbm;
            hbz = rbz / rbm;
        }

        double axbx = hay * hbz - haz * hby;
        double axby = haz * hbx - hax * hbz;
        double axbz = hax * hby - hay * hbx;

        double bxcy = -hbz;
        double bxcz = hby;

        double cxay = haz;
        double cxaz = -hay;

        double abden = 1.0 + hax * hbx + hay * hby + haz * hbz;
        double abfac = 0.0;
        if (abden > tol) {
            abfac = (rar + rbr) / abden;
        }

        double bcden = 1.0 - hbx;
        double bcfac = 0.0;
        if (bcden > tol) {
            bcfac = rbr / bcden;
        }

        double caden = 1.0 - hax;
        double cafac = 0.0;
        if (caden > tol) {
            cafac = rar / caden;
        }

        vx = (axbx * abfac)/fourpi;
        vy = (axby * abfac + bxcy * bcfac + cxay * cafac)/fourpi;
        vz = (axbz * abfac + bxcz * bcfac + cxaz * cafac)/fourpi;

        ''', 'cupy_cwdv')


    def cupy_cwdv(pnts: Vector, grda: Vector, grdb: Vector,
                *, tol: float = 1e-12, time: bool = False,
                betm: float = 1.0) -> Vector:

        if time:
            start1 = perf_counter()
        px, py, pz = asarray(pnts.x), asarray(pnts.y), asarray(pnts.z)
        ax, ay, az = asarray(grda.x), asarray(grda.y), asarray(grda.z)
        bx, by, bz = asarray(grdb.x), asarray(grdb.y), asarray(grdb.z)
        if time:
            finish1 = perf_counter()
            elapsed1 = finish1 - start1
            print(f'Cupy in time elapsed is {elapsed1:.6f} seconds.')

        if time:
            start2 = perf_counter()
        vx, vy, vz = cupy_cwdv_kernel(px, py, pz, ax, ay, az, bx, by, bz,
                                      tol, betm)
        if time:
            finish2 = perf_counter()
            elapsed2 = finish2 - start2
            print(f'Cupy execution time elapsed is {elapsed2:.6f} seconds.')

        if time:
            start3 = perf_counter()
        vel = Vector(vx.get(), vy.get(), vz.get())
        if time:
            finish3 = perf_counter()
            elapsed3 = finish3 - start3
            print(f'Cupy out time elapsed is {elapsed3:.6f} seconds.')

        return vel

    # Example usage of the cupy_cwdv function
    pnts = Vector.zeros((1, ))
    grda = Vector(-1.2, 1.0, 0.0)
    grdb = Vector(-0.8, -1.0, 0.0)

    _ = cupy_cwdv(pnts, grda, grdb)

except ImportError:
    pass
