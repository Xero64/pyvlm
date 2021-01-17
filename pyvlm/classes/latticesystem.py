from math import pi
from numpy.matlib import zeros, multiply, divide, fill_diagonal, seterr, logical_and
from pygeom.geom3d import Point
from pygeom.matrix3d import zero_matrix_vector, MatrixVector
from pygeom.matrix3d import elementwise_dot_product, elementwise_cross_product
from pygeom.matrix3d import elementwise_multiply, elementwise_divide
from py2md.classes import MDTable

seterr(divide='ignore', invalid='ignore')

fourPi = 4*pi

class LatticeSystem(object):
    source = None # System Source
    name = None # System Name
    srfcs = None # System Surfaces
    strps = None # System Strips
    pnls = None # System Panels
    bref = None # Reference Span
    cref = None # Reference Chord
    sref = None # Reference Area
    rref = None # Reference Point
    ctrls = None # System Controls
    nump = None # Number of Panels
    nums = None # Number of Strips
    masses = None # Store Mass Options
    _ra = None # Horseshoe Vortex Point A
    _rb = None # Horseshoe Vortex Point B
    _rc = None # Panel Control Point
    _rg = None # Panel Induced Point
    _ungam = None # System Unit Solution
    _avg = None # Induced Point Velocity Matrix
    _afg = None # Induced Point Force Matrix
    _avc = None # Control Point Velocity Matrix
    _aic = None # Influence Coefficient Matrix
    _afs = None
    _adc = None
    _ada = None
    _ava = None
    _avb = None
    _bvg = None
    _bdg = None
    _blg = None
    _byg = None
    _bmg = None
    _bda = None
    _ar = None # Aspect Ratio
    _cdo = None
    _cdo_ff = None
    results = None
    def __init__(self, name: str, srfcs: list,
                       bref: float, cref: float, sref: float,
                       rref: Point):
        self.name = name
        self.srfcs = srfcs
        self.bref = bref
        self.cref = cref
        self.sref = sref
        self.rref = rref
        self.results = {}
    def mesh(self):
        lsid = 0
        lpid = 0
        for srfc in self.srfcs:
            lsid, lpid = srfc.mesh(lsid, lpid)
        self.strps = []
        self.pnls = []
        self.ctrls = {}
        for srfc in self.srfcs:
            for strp in srfc.strps:
                self.strps.append(strp)
            pnls = srfc.return_panels()
            for pnl in pnls:
                self.pnls.append(pnl)
        self.nump = len(self.pnls)
        self.nums = len(self.strps)
        ind = 2
        for srfc in self.srfcs:
            for sht in srfc.shts:
                for control in sht.ctrls:
                    if control not in self.ctrls:
                        self.ctrls[control] = (ind, ind+1, ind+2, ind+3)
                        ind += 4
    def reset(self):
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None
    @property
    def ra(self):
        if self._ra is None:
            self._ra = zero_matrix_vector((1, self.nump))
            for pnl in self.pnls:
                self._ra[0, pnl.lpid] = pnl.pnta
        return self._ra
    @property
    def rb(self):
        if self._rb is None:
            self._rb = zero_matrix_vector((1, self.nump))
            for pnl in self.pnls:
                self._rb[0, pnl.lpid] = pnl.pntb
        return self._rb
    @property
    def rc(self):
        if self._rc is None:
            self._rc = zero_matrix_vector((self.nump, 1))
            for pnl in self.pnls:
                self._rc[pnl.lpid, 0] = pnl.pntc
        return self._rc
    @property
    def rg(self):
        if self._rg is None:
            self._rg = zero_matrix_vector((self.nump, 1))
            for pnl in self.pnls:
                self._rg[pnl.lpid, 0] = pnl.pnti
        return self._rg
    def avc(self, mach: float):
        if self._avc is None:
            self._avc = {}
        if mach not in self._avc:
            beta = (1-mach**2)**0.5
            veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rc, beta)
            self._avc[mach] = (veli+vela-velb)/fourPi
        return self._avc[mach]
    def aic(self, mach: float):
        if self._aic is None:
            self._aic = {}
        if mach not in self._aic:
            avc = self.avc(mach)
            aic = zeros(avc.shape)
            for pnl in self.pnls:
                aic[pnl.lpid, :] = avc[pnl.lpid, :]*pnl.nrml
            self._aic[mach] = aic
        return self._aic[mach]
    @property
    def afs(self):
        if self._afs is None:
            num = len(self.pnls)
            numc = len(self.ctrls)
            self._afs = zero_matrix_vector((num, 2+4*numc))
            for pnl in self.pnls:
                lpid = pnl.lpid
                rrel = pnl.pntc-self.rref
                self._afs[lpid, 0] = pnl.nrml
                self._afs[lpid, 1] = -rrel**pnl.nrml
            for srfc in self.srfcs:
                for sht in srfc.shts:
                    for control in sht.ctrls:
                        ctrl = sht.ctrls[control]
                        ctup = self.ctrls[control]
                        for pnl in ctrl.pnls:
                            lpid = pnl.lpid
                            rrel = pnl.pntc-self.rref
                            dndlp = pnl.dndl(ctrl.posgain, ctrl.uhvec)
                            self._afs[lpid, ctup[0]] = dndlp
                            self._afs[lpid, ctup[1]] = rrel**dndlp
                            dndln = pnl.dndl(ctrl.neggain, ctrl.uhvec)
                            self._afs[lpid, ctup[2]] = dndln
                            self._afs[lpid, ctup[3]] = rrel**dndln
        return self._afs
    def ungam(self, mach: float):
        if self._ungam is None:
            self._ungam = {}
        if mach not in self._ungam:
            from pygeom.matrix3d import solve_matrix_vector
            aic = self.aic(mach)
            self._ungam[mach] = -solve_matrix_vector(aic, self.afs)
        return self._ungam[mach]
    def avg(self, mach: float):
        if self._avg is None:
            self._avg = {}
        if mach not in self._avg:
            beta = (1-mach**2)**0.5
            veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rg, beta)
            fill_diagonal(veli.x, 0.0)
            fill_diagonal(veli.y, 0.0)
            fill_diagonal(veli.z, 0.0)
            self._avg[mach] = (veli+vela-velb)/fourPi
        return self._avg[mach]
    def afg(self, mach: float):
        if self._afg is None:
            self._afg = {}
        if mach not in self._afg:
            avg = self.avg(mach)
            afg = zero_matrix_vector(avg.shape)
            for pnl in self.pnls:
                if not pnl.noload:
                    afg[pnl.lpid, :] = avg[pnl.lpid, :]**pnl.leni
            self._afg[mach] = afg
        return self._afg[mach]
    def adc(self, mach: float):
        if self._adc is None:
            self._adc = {}
        if mach not in self._adc:
            avc = self.avc(mach)
            adc = zeros(avc.shape)
            for pnl in self.pnls:
                adc[pnl.lpid, :] = avc[pnl.lpid, :]*pnl.tang
            self._adc[mach] = adc
        return self._adc[mach]
    @property
    def ada(self):
        if self._ada is None:
            num = len(self.pnls)
            self._ada = zeros((num, 1))
            for pnl in self.pnls:
                self._ada[pnl.lpid, 0] = pnl.cdoarea
        return self._ada
    @property
    def cdo(self):
        if self._cdo is None:
            dragarea = 0.0
            for pnl in self.pnls:
                dragarea += self.ada[pnl.lpid, 0]
            self._cdo = dragarea/self.sref
        return self._cdo
    @property
    def bvg(self):
        if self._bvg is None:
            num = len(self.strps)
            self._bvg = zeros((num, num))
            for strpi in self.strps:
                for strpj in self.strps:
                    self._bvg[strpi.lsid, strpj.lsid] = strpj.trefftz_velocity(strpi.pnti)*strpi.nrmt
        return self._bvg
    @property
    def bdg(self):
        if self._bdg is None:
            num = len(self.strps)
            self._bdg = zeros((num, num))
            for strpi in self.strps:
                for strpj in self.strps:
                    self._bdg[strpi.lsid, strpj.lsid] = strpi.trefftz_drag(self.bvg[strpi.lsid, strpj.lsid])
        return self._bdg
    @property
    def blg(self):
        if self._blg is None:
            num = len(self.strps)
            self._blg = zeros((num, 1))
            for strp in self.strps:
                self._blg[strp.lsid, 0] = strp.trefftz_lift()
        return self._blg
    @property
    def byg(self):
        if self._byg is None:
            num = len(self.strps)
            self._byg = zeros((num, 1))
            for strp in self.strps:
                self._byg[strp.lsid, 0] = strp.trefftz_yfrc()
        return self._byg
    @property
    def bmg(self):
        if self._bmg is None:
            num = len(self.strps)
            self._bmg = zeros((num, 1))
            for strp in self.strps:
                self._bmg[strp.lsid, 0] += strp.pnti.y*self.blg[strp.lsid, 0]
                self._bmg[strp.lsid, 0] -= strp.pnti.z*self.byg[strp.lsid, 0]
        return self._bmg
    @property
    def bda(self):
        if self._bda is None:
            num = len(self.strps)
            self._bda = zeros((num, 1))
            for strp in self.strps:
                self._bda[strp.lsid, 0] = strp.cdoarea
        return self._bda
    @property
    def cdo_ff(self):
        if self._cdo_ff is None:
            dragarea = 0.0
            for strp in self.strps:
                dragarea += self.bda[strp.lsid, 0]
            self._cdo_ff = dragarea/self.sref
        return self._cdo_ff
    @property
    def ar(self):
        if self._ar is None:
            self._ar = self.bref**2/self.sref
        return self._ar
    @property
    def lstrpi(self):
        sgrp = []
        for srfc in self.srfcs:
            sgrp += srfc.sgrp[0]
        return sgrp
    @property
    def mstrpi(self):
        sgrp = []
        for srfc in self.srfcs:
            sgrp += srfc.sgrp[1]
        return sgrp
    def set_strip_alpha(self, alpha: list):
        for strp in self.strps:
            strp.set_angle(alpha[strp.lsid])
        self._aic = None
        self._afs = None
        self._ungam = None
    @property
    def strip_geometry(self):
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('xpos', '.5f')
        table.add_column('ypos', '.5f')
        table.add_column('Zle', '.5f')
        table.add_column('Chord', '.4f')
        table.add_column('Width', '.5f')
        table.add_column('Area', '.6f')
        table.add_column('Dihed', '.4f')
        table.add_column('Incid', '.4f')
        for strp in self.strps:
            j = strp.lsid
            xpos = strp.pnti.x
            ypos = strp.pnti.y
            zpos = strp.pnti.z
            chord = strp.chord
            width = strp.dst
            area = strp.area
            dihed = strp.dihedral
            angle = strp.angle
            table.add_row([j, xpos, ypos, zpos, chord, width, area, dihed, angle])
        return table
    @property
    def panel_geometry(self):
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('X', '.5f')
        table.add_column('Y', '.5f')
        table.add_column('Z', '.5f')
        table.add_column('DX', '.5f')
        table.add_column('nx', '.5f')
        table.add_column('ny', '.5f')
        table.add_column('nz', '.5f')
        for pnl in self.pnls:
            j = pnl.lpid
            x = pnl.pntg.x
            y = pnl.pntg.y
            z = pnl.pntg.z
            dx = pnl.crd
            nx = pnl.nrml.x
            ny = pnl.nrml.y
            nz = pnl.nrml.z
            table.add_row([j, x, y, z, dx, nx, ny, nz])
        return table
    def copy_from_source(self):
        lsys = latticesystem_from_json(self.source)
        for attr in self.__dict__:
            if attr[0] == '_':
                lsys.__dict__[attr] = self.__dict__[attr].copy()
        return lsys
    def velocity_matrix(self, rc: MatrixVector):
        veli, vela, velb = velocity_matrix(self.ra, self.rb, rc)
        return (veli+vela-velb)/fourPi
    def __repr__(self):
        return '<LatticeSystem: {:s}>'.format(self.name)
    def __str__(self):
        outstr = '# Lattice System '+self.name+'\n'
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('Sref', 'g', data=[self.sref])
        table.add_column('cref', 'g', data=[self.cref])
        table.add_column('bref', 'g', data=[self.bref])
        table.add_column('xref', '.3f', data=[self.rref.x])
        table.add_column('yref', '.3f', data=[self.rref.y])
        table.add_column('zref', '.3f', data=[self.rref.z])
        outstr += table._repr_markdown_()
        table = MDTable()
        if self.strps is not None:
            table.add_column('# Strips', 'd', data=[len(self.strps)])
        if self.pnls is not None:
            table.add_column('# Panels', 'd', data=[len(self.pnls)])
        if self.ctrls is not None:
            table.add_column('# Controls', 'd', data=[len(self.ctrls)])
        if len(table.columns) > 0:
            outstr += table._repr_markdown_()
        return outstr
    def _repr_markdown_(self):
        return self.__str__()

def latticesystem_from_json(jsonfilepath: str, mesh: bool=True):
    from json import load

    with open(jsonfilepath, 'rt') as jsonfile:
        sysdct = load(jsonfile)

    sysdct['source'] = jsonfilepath

    sys = latticesystem_from_dict(sysdct)

    if mesh and sys.pnls is None:
        sys.mesh()

    return sys

def latticesystem_from_dict(sysdct: dict):
    from .latticesurface import latticesurface_from_json
    from .latticeresult import latticeresult_from_json
    from .latticetrim import latticetrim_from_json
    from pyvlm.tools import masses_from_json, masses_from_data
    from os.path import dirname, join, exists

    jsonfilepath = sysdct['source']

    path = dirname(jsonfilepath)

    for surfdata in sysdct['surfaces']:
        for sectdata in surfdata['sections']:
            if 'airfoil' in sectdata:
                airfoil = sectdata['airfoil']
                if airfoil[-4:] == '.dat':
                    airfoil = join(path, airfoil)
                    if not exists(airfoil):
                        print(f'Airfoil {airfoil} does not exist.')
                        del sectdata['airfoil']
                    else:
                        sectdata['airfoil'] = airfoil

    name = sysdct['name']
    sfcs = []
    for surfdata in sysdct['surfaces']:
        sfc = latticesurface_from_json(surfdata)
        sfcs.append(sfc)
    bref = sysdct['bref']
    cref = sysdct['cref']
    sref = sysdct['sref']
    xref = sysdct['xref']
    yref = sysdct['yref']
    zref = sysdct['zref']
    rref = Point(xref, yref, zref)
    lsys = LatticeSystem(name, sfcs, bref, cref, sref, rref)

    masses = {}
    if 'masses' in sysdct:
        if isinstance(sysdct['masses'], list):
            masses = masses_from_data(sysdct['masses'])
        elif isinstance(sysdct['masses'], str):
            if sysdct['masses'][-5:] == '.json':
                massfilename = sysdct['masses']
                massfilepath = join(path, massfilename)
            masses = masses_from_json(massfilepath)
    lsys.masses = masses

    if 'cases' in sysdct and sysdct:
        lsys.mesh()
        for i in range(len(sysdct['cases'])):
            resdata = sysdct['cases'][i]
            if 'trim' in resdata:
                latticetrim_from_json(lsys, resdata)
            else:
                latticeresult_from_json(lsys, resdata)

    lsys.source = jsonfilepath

    return lsys

def velocity_matrix(ra: MatrixVector, rb: MatrixVector, rc: MatrixVector,
                    betm: float=1.0, tol: float=1e-12):

    if ra.shape != rb.shape:
        return ValueError()

    numi = rc.shape[0]
    numj = ra.shape[1]

    ra = ra.repeat(numi, axis=0)
    rb = rb.repeat(numi, axis=0)
    rc = rc.repeat(numj, axis=1)

    a = rc-ra
    b = rc-rb

    a.x = a.x/betm
    b.x = b.x/betm

    am = a.return_magnitude()
    bm = b.return_magnitude()

    # Velocity from Bound Vortex
    adb = elementwise_dot_product(a, b)
    abm = multiply(am, bm)
    dm = multiply(abm, abm+adb)
    axb = elementwise_cross_product(a, b)
    axbm = axb.return_magnitude()
    chki = (axbm == 0.0)
    chki = logical_and(axbm >= -tol, axbm <= tol)
    veli = elementwise_multiply(axb, divide(am+bm, dm))
    veli.x[chki] = 0.0
    veli.y[chki] = 0.0
    veli.z[chki] = 0.0

    # Velocity from Trailing Vortex A
    axx = MatrixVector(zeros(a.shape), a.z, -a.y)
    axxm = axx.return_magnitude()
    chka = (axxm == 0.0)
    vela = elementwise_divide(axx, multiply(am, am-a.x))
    vela.x[chka] = 0.0
    vela.y[chka] = 0.0
    vela.z[chka] = 0.0

    # Velocity from Trailing Vortex B
    bxx = MatrixVector(zeros(b.shape), b.z, -b.y)
    bxxm = bxx.return_magnitude()
    chkb = (bxxm == 0.0)
    velb = elementwise_divide(bxx, multiply(bm, bm-b.x))
    velb.x[chkb] = 0.0
    velb.y[chkb] = 0.0
    velb.z[chkb] = 0.0

    return veli, vela, velb
