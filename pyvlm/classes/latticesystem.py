from typing import TYPE_CHECKING, Dict, List, Tuple, Union

from numpy import absolute, divide, multiply, pi, reciprocal, sqrt, zeros
from py2md.classes import MDTable
from pygeom.array3d import ArrayVector, solve_arrayvector, zero_arrayvector
from pygeom.geom3d import Vector

if TYPE_CHECKING:
    from numpy import float64
    from numpy.typing import NDArray

    from ..tools.mass import Mass, MassCollection
    from .latticepanel import LatticePanel
    from .latticeresult import LatticeResult
    from .latticestrip import LatticeStrip
    from .latticesurface import LatticeSurface

    MassLike = Union[Mass, MassCollection]

FOURPI = 4*pi


class LatticeSystem():
    source: str = None # System Source
    name: str = None # System Name
    srfcs: List['LatticeSurface'] = None # System Surfaces
    strps: List['LatticeStrip'] = None # System Strips
    pnls: List['LatticePanel'] = None # System Panels
    bref: float = None # Reference Span
    cref: float = None # Reference Chord
    sref: float = None # Reference Area
    rref: Vector = None # Reference Vector
    ctrls: Dict[str, Tuple[int, int, int, int]] = None # System Controls
    nump: int = None # Number of Panels
    nums: int = None # Number of Strips
    masses: Dict[str, 'MassLike'] = None # Store Mass Options
    _ra: ArrayVector = None # Horseshoe Vortex Vector A
    _rb: ArrayVector = None # Horseshoe Vortex Vector B
    _rc: ArrayVector = None # Panel Control Vector
    _rg: ArrayVector = None # Panel Induced Vector
    _ungam: ArrayVector = None # System Unit Solution
    _avg: ArrayVector = None # Induced Vector Velocity Matrix
    _afg: ArrayVector = None # Induced Vector Force Matrix
    _avc: ArrayVector = None # Control Vector Velocity Matrix
    _aic: 'NDArray[float64]' = None # Influence Coefficient Matrix
    _afs: ArrayVector = None
    _ada: 'NDArray[float64]' = None
    _bvg: ArrayVector = None
    _bdg: 'NDArray[float64]' = None
    _blg: 'NDArray[float64]' = None
    _byg: 'NDArray[float64]' = None
    _bmg: 'NDArray[float64]' = None
    _bda: 'NDArray[float64]' = None
    _ar: float = None # Aspect Ratio
    _cdo: float = None
    _cdo_ff: float = None
    results: Dict[str, 'LatticeResult'] = None

    def __init__(self, name: str, srfcs: list,
                       bref: float, cref: float, sref: float,
                       rref: Vector):
        self.name = name
        self.srfcs = srfcs
        self.bref = bref
        self.cref = cref
        self.sref = sref
        self.rref = rref
        self.results = {}

    def mesh(self) -> None:
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

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr[0] == '_':
                self.__dict__[attr] = None

    @property
    def ra(self) -> ArrayVector:
        if self._ra is None:
            self._ra = zero_arrayvector(self.nump, dtype=float)
            for pnl in self.pnls:
                self._ra[pnl.lpid] = pnl.pnta
        return self._ra

    @property
    def rb(self) -> ArrayVector:
        if self._rb is None:
            self._rb = zero_arrayvector(self.nump, dtype=float)
            for pnl in self.pnls:
                self._rb[pnl.lpid] = pnl.pntb
        return self._rb

    @property
    def rc(self) -> ArrayVector:
        if self._rc is None:
            self._rc = zero_arrayvector(self.nump, dtype=float)
            for pnl in self.pnls:
                self._rc[pnl.lpid] = pnl.pntc
        return self._rc

    @property
    def rg(self) -> ArrayVector:
        if self._rg is None:
            self._rg = zero_arrayvector(self.nump, dtype=float)
            for pnl in self.pnls:
                self._rg[pnl.lpid] = pnl.pnti
        return self._rg

    def avc(self, mach: float) -> ArrayVector:
        if self._avc is None:
            self._avc = {}
        if mach not in self._avc:
            beta = sqrt(1.0 - mach**2)
            veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rc, beta)
            self._avc[mach] = (veli + vela - velb)/FOURPI
        return self._avc[mach]

    def aic(self, mach: float) -> 'NDArray[float64]':
        if self._aic is None:
            self._aic = {}
        if mach not in self._aic:
            avc = self.avc(mach)
            aic = zeros(avc.shape, dtype=float)
            for pnl in self.pnls:
                aic[pnl.lpid, :] = avc[pnl.lpid, :].dot(pnl.nrml)
            self._aic[mach] = aic
        return self._aic[mach]

    @property
    def afs(self) -> ArrayVector:
        if self._afs is None:
            num = len(self.pnls)
            numc = len(self.ctrls)
            self._afs = zero_arrayvector((num, 2+4*numc), dtype=float)
            for pnl in self.pnls:
                lpid = pnl.lpid
                rrel = pnl.pntc-self.rref
                self._afs[lpid, 0] = pnl.nrml
                self._afs[lpid, 1] = -rrel.cross(pnl.nrml)
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
                            self._afs[lpid, ctup[1]] = rrel.cross(dndlp)
                            dndln = pnl.dndl(ctrl.neggain, ctrl.uhvec)
                            self._afs[lpid, ctup[2]] = dndln
                            self._afs[lpid, ctup[3]] = rrel.cross(dndln)
        return self._afs

    def ungam(self, mach: float) -> ArrayVector:
        if self._ungam is None:
            self._ungam = {}
        if mach not in self._ungam:
            aic = self.aic(mach)
            self._ungam[mach] = -solve_arrayvector(aic, self.afs)
        return self._ungam[mach]

    def avg(self, mach: float) -> ArrayVector:
        if self._avg is None:
            self._avg = {}
        if mach not in self._avg:
            beta = (1 - mach**2)**0.5
            veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rg, beta)
            self._avg[mach] = (veli + vela - velb)/FOURPI
        return self._avg[mach]

    def afg(self, mach: float) -> ArrayVector:
        if self._afg is None:
            self._afg = {}
        if mach not in self._afg:
            avg = self.avg(mach)
            afg = zero_arrayvector(avg.shape, dtype=float)
            for pnl in self.pnls:
                if not pnl.noload:
                    afg[pnl.lpid, :] = avg[pnl.lpid, :].cross(pnl.leni)
            self._afg[mach] = afg
        return self._afg[mach]

    @property
    def ada(self) -> 'NDArray[float64]':
        if self._ada is None:
            num = len(self.pnls)
            self._ada = zeros(num, dtype=float)
            for pnl in self.pnls:
                self._ada[pnl.lpid] = pnl.cdoarea
        return self._ada

    @property
    def cdo(self) -> 'NDArray[float64]':
        if self._cdo is None:
            dragarea = 0.0
            for pnl in self.pnls:
                dragarea += self.ada[pnl.lpid]
            self._cdo = dragarea/self.sref
        return self._cdo

    @property
    def bvg(self) -> 'NDArray[float64]':
        if self._bvg is None:
            num = len(self.strps)
            self._bvg = zeros((num, num), dtype=float)
            for strpi in self.strps:
                for strpj in self.strps:
                    self._bvg[strpi.lsid, strpj.lsid] = strpj.trefftz_velocity(strpi.pnti).dot(strpi.nrmt)
        return self._bvg

    @property
    def bdg(self) -> 'NDArray[float64]':
        if self._bdg is None:
            num = len(self.strps)
            self._bdg = zeros((num, num), dtype=float)
            for strpi in self.strps:
                for strpj in self.strps:
                    self._bdg[strpi.lsid, strpj.lsid] = strpi.trefftz_drag(self.bvg[strpi.lsid, strpj.lsid])
        return self._bdg

    @property
    def blg(self) -> 'NDArray[float64]':
        if self._blg is None:
            num = len(self.strps)
            self._blg = zeros(num, dtype=float)
            for strp in self.strps:
                self._blg[strp.lsid] = strp.trefftz_lift()
        return self._blg

    @property
    def byg(self) -> 'NDArray[float64]':
        if self._byg is None:
            num = len(self.strps)
            self._byg = zeros(num, dtype=float)
            for strp in self.strps:
                self._byg[strp.lsid] = strp.trefftz_yfrc()
        return self._byg

    @property
    def bmg(self) -> 'NDArray[float64]':
        if self._bmg is None:
            num = len(self.strps)
            self._bmg = zeros(num, dtype=float)
            for strp in self.strps:
                self._bmg[strp.lsid] += strp.pnti.y*self.blg[strp.lsid]
                self._bmg[strp.lsid] -= strp.pnti.z*self.byg[strp.lsid]
        return self._bmg

    @property
    def bda(self) -> 'NDArray[float64]':
        if self._bda is None:
            num = len(self.strps)
            self._bda = zeros(num, dtype=float)
            for strp in self.strps:
                self._bda[strp.lsid] = strp.cdoarea
        return self._bda

    @property
    def cdo_ff(self) -> float:
        if self._cdo_ff is None:
            dragarea = 0.0
            for strp in self.strps:
                dragarea += self.bda[strp.lsid]
            self._cdo_ff = dragarea/self.sref
        return self._cdo_ff

    @property
    def ar(self) -> float:
        if self._ar is None:
            self._ar = self.bref**2/self.sref
        return self._ar

    @property
    def lstrpi(self) -> List['LatticeStrip']:
        sgrp = []
        for srfc in self.srfcs:
            sgrp += srfc.sgrp[0]
        return sgrp

    @property
    def mstrpi(self) -> List['LatticeStrip']:
        sgrp = []
        for srfc in self.srfcs:
            sgrp += srfc.sgrp[1]
        return sgrp

    def set_strip_alpha(self, alpha: List[float]) -> None:
        for strp in self.strps:
            strp.set_twist(alpha[strp.lsid])
        self._aic = None
        self._afs = None
        self._ungam = None

    @property
    def strip_geometry(self) -> MDTable:
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
            twist = strp.twist
            table.add_row([j, xpos, ypos, zpos, chord, width, area, dihed, twist])
        return table

    @property
    def panel_geometry(self) -> MDTable:
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

    def copy_from_source(self) -> 'LatticeSystem':
        lsys = latticesystem_from_json(self.source)
        for attr in self.__dict__:
            if attr[0] == '_':
                lsys.__dict__[attr] = self.__dict__[attr].copy()
        return lsys

    def velocity_matrix(self, rc: ArrayVector) -> ArrayVector:
        veli, vela, velb = velocity_matrix(self.ra, self.rb, rc)
        return (veli + vela - velb)/FOURPI

    def __repr__(self) -> str:
        return '<LatticeSystem: {:s}>'.format(self.name)

    def __str__(self) -> str:
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

    def _repr_markdown_(self) -> str:
        return self.__str__()


def latticesystem_from_json(jsonfilepath: str, mesh: bool=True) -> LatticeSystem:
    from json import load

    with open(jsonfilepath, 'rt') as jsonfile:
        sysdct = load(jsonfile)

    sysdct['source'] = jsonfilepath

    sys = latticesystem_from_dict(sysdct)

    if mesh and sys.pnls is None:
        sys.mesh()

    return sys

def latticesystem_from_dict(sysdct: dict) -> LatticeSystem:
    from os.path import dirname, exists, join

    from pyvlm.tools import masses_from_data, masses_from_json

    from .latticeresult import latticeresult_from_json
    from .latticesurface import latticesurface_from_json
    from .latticetrim import latticetrim_from_json

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
    rref = Vector(xref, yref, zref)
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

def velocity_matrix(ra: ArrayVector, rb: ArrayVector, rc: ArrayVector,
                    betm: float=1.0, tol: float=1e-12) -> ArrayVector:

    if ra.size != rb.size:
        raise ValueError('ra and rb sizes do not match.')

    numi = rc.size
    numj = ra.size

    ra = ra.reshape((1, -1)).repeat(numi, axis=0)
    rb = rb.reshape((1, -1)).repeat(numi, axis=0)
    rc = rc.reshape((-1, 1)).repeat(numj, axis=1)

    a = rc - ra
    b = rc - rb

    a.x = a.x/betm
    b.x = b.x/betm

    am = a.return_magnitude()
    bm = b.return_magnitude()

    # Velocity from Bound Vortex
    adb = a.dot(b)
    abm = multiply(am, bm)

    # from numpy import argwhere

    # dmi_min = dmi.min()

    # ind_min = argwhere(dmi == dmi_min)

    # print(f'dmi_min = {dmi_min}')
    # print(f'ind_min = {ind_min}')

    axb = a.cross(b)
    dmi = multiply(abm, abm + adb)
    chki = absolute(dmi) > tol
    faci = zeros(dmi.shape)
    divide(am + bm, dmi, where=chki, out=faci)
    veli = axb*faci
    # veli.x[chki] = 0.0
    # veli.y[chki] = 0.0
    # veli.z[chki] = 0.0

    # Velocity from Trailing Vortex A
    axx = ArrayVector(zeros(a.shape, dtype=float), a.z, -a.y)
    dma = multiply(am, am - a.x)
    chka = absolute(dma) > tol
    faca = zeros(dma.shape)
    reciprocal(dma, where=chka, out=faca)
    vela = axx*faca
    # vela.x[chka] = 0.0
    # vela.y[chka] = 0.0
    # vela.z[chka] = 0.0

    # Velocity from Trailing Vortex B
    bxx = ArrayVector(zeros(b.shape, dtype=float), b.z, -b.y)
    dmb = multiply(bm, bm - b.x)
    chkb = absolute(dmb) > tol
    facb = zeros(dmb.shape)
    reciprocal(dmb, where=chkb, out=facb)
    velb = bxx*facb
    # velb.x[chkb] = 0.0
    # velb.y[chkb] = 0.0
    # velb.z[chkb] = 0.0

    return veli, vela, velb
