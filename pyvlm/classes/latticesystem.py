from time import perf_counter
from math import pi, atan2
from numpy.matlib import zeros, matrix, multiply, divide, fill_diagonal, seterr
from pygeom.geom3d import Point, Vector, ihat, jhat, zero_vector
from pygeom.matrix3d import zero_matrix_vector, MatrixVector
from pygeom.matrix3d import elementwise_dot_product, elementwise_cross_product
from pygeom.matrix3d import elementwise_multiply, elementwise_divide

seterr(divide='ignore', invalid='ignore')

class LatticeSystem(object):
    source = None # System Source
    name = None # System Name
    srfcs = None # System Surfaces
    strps = None # System Strips
    mstrp = None # Mirrored Strips
    pnls = None # System Panels
    mpnl = None # Mirrored Panels
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
    # _rt = None
    # _rtva = None
    # _rtvb = None
    _gam = None # System Unit Solution
    _avg = None # Induced Point Velocity Matrix
    _afg = None # Induced Point Force Matrix
    # _ava = None
    # _afa = None
    # _avb = None
    # _afb = None
    # _avt = None
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
        self.mesh()
        self.inherit()
        # self.build()
    def mesh(self):
        lsid = 0
        lpid = 0
        for srfc in self.srfcs:
            lsid, lpid = srfc.mesh(lsid, lpid)
    def inherit(self):
        strpdct = {}
        pnldct = {}
        istrp = 0
        ipnl = 0
        self.strps = []
        self.mstrp = []
        self.pnls = []
        self.mpnl = []
        self.ctrls = {}
        for srfc in self.srfcs:
            for strp in srfc.strps:
                self.strps.append(strp)
                if strp.msid is None:
                    strpind = istrp
                    strpdct[strp.lsid] = istrp
                else:
                    strpind = strpdct[strp.msid]
                self.mstrp.append(strpind)
                istrp += 1
            pnls = srfc.return_panels()
            for pnl in pnls:
                self.pnls.append(pnl)
                if pnl.mpid is None:
                    pnlind = ipnl
                    pnldct[pnl.lpid] = ipnl
                else:
                    pnlind = pnldct[pnl.mpid]
                self.mpnl.append(pnlind)
                ipnl += 1
        self.nump = len(self.pnls)
        self.nums = len(self.strps)
        ind = 2
        for srfc in self.srfcs:
            for sht in srfc.shts:
                for control in sht.ctrls:
                    if control not in self.ctrls:
                        self.ctrls[control] = (ind, ind+1, ind+2, ind+3)
                        ind += 4
    def build(self):
        self.ra
        self.rb
        self.rc
        self.rg
        self.avc
        self.aic
        self.afs
        self.gam
        self.avg
        self.afg
        # self.ava
        # self.afa
        # self.avb
        # self.afb
        self.bvg
        self.bdg
        self.blg
        self.byg
        self.bmg
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
    # @property
    # def rt(self):
    #     if self._rt is None:
    #         num = len(self.strps)
    #         self._rt = zero_matrix_vector((num, 1))
    #         for strpi in self.strps:
    #             i = strpi.lsid
    #             self._rt[i, :] = strpi.pnte
    #     return self._rt
    # @property
    # def rtva(self):
    #     if self._rtva is None:
    #         num = len(self.pnls)
    #         self._rtva = zero_matrix_vector((num, num))
    #         for pnli in self.pnls:
    #             i = pnli.lpid
    #             self._rtva[i, :] = pnli.ptva
    #     return self._rtva
    # @property
    # def rtvb(self):
    #     if self._rtvb is None:
    #         num = len(self.pnls)
    #         self._rtvb = zero_matrix_vector((num, num))
    #         for pnli in self.pnls:
    #             i = pnli.lpid
    #             self._rtvb[i, :] = pnli.ptvb
    #     return self._rtvb
    @property
    def avc(self):
        if self._avc is None:
            veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rc)
            self._avc = (veli+vela-velb)/(4*pi)
        return self._avc
    @property
    def aic(self):
        if self._aic is None:
            self._aic = zeros((self.nump, self.nump))
            for pnl in self.pnls:
                self._aic[pnl.lpid, :] = self.avc[pnl.lpid, :]*pnl.nrml
        return self._aic
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
    @property
    def gam(self):
        if self._gam is None:
            from pygeom.matrix3d import solve_matrix_vector
            self._gam = -solve_matrix_vector(self.aic, self.afs)
        return self._gam
    @property
    def avg(self):
        if self._avg is None:
            veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rg)
            fill_diagonal(veli.x, 0.0)
            fill_diagonal(veli.y, 0.0)
            fill_diagonal(veli.z, 0.0)
            self._avg = (veli+vela-velb)/(4*pi)
        return self._avg
    # @property
    # def avt(self):
    #     if self._avt is None:
    #         veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rt)
    #         self._avt = (veli+vela-velb)/(4*pi)
    #     return self._avt
    @property
    def afg(self):
        if self._afg is None:
            num = len(self.pnls)
            self._afg = zero_matrix_vector((num, num))
            for pnl in self.pnls:
                if not pnl.noload:
                    i = pnl.lpid
                    self._afg[i, :] = self.avg[i, :]**pnl.leni
        return self._afg
    # @property
    # def ava(self):
    #     if self._ava is None:
    #         veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rtva)
    #         fill_diagonal(vela.x, 0.0)
    #         fill_diagonal(vela.y, 0.0)
    #         fill_diagonal(vela.z, 0.0)
    #         self._ava = (veli+vela-velb)/(4*pi)
    #     return self._ava
    # @property
    # def afa(self):
    #     if self._afa is None:
    #         num = len(self.pnls)
    #         self._afa = zero_matrix_vector((num, num))
    #         for pnl in self.pnls:
    #             if not pnl.noload:
    #                 i = pnl.lpid
    #                 self._afa[i, :] = self.ava[i, :]**pnl.ltva
    #     return self._afa
    # @property
    # def avb(self):
    #     if self._avb is None:
    #         veli, vela, velb = velocity_matrix(self.ra, self.rb, self.rtvb)
    #         fill_diagonal(velb.x, 0.0)
    #         fill_diagonal(velb.y, 0.0)
    #         fill_diagonal(velb.z, 0.0)
    #         self._avb = (veli+vela-velb)/(4*pi)
    #     return self._avb
    # @property
    # def afb(self):
    #     if self._afb is None:
    #         num = len(self.pnls)
    #         self._afb = zero_matrix_vector((num, num))
    #         for pnl in self.pnls:
    #             if not pnl.noload:
    #                 i = pnl.lpid
    #                 self._afb[i, :] = self.avb[i, :]**pnl.ltvb
    #     return self._afb
    @property
    def adc(self):
        if self._adc is None:
            num = len(self.pnls)
            self._adc = zeros((num, num))
            for pnl in self.pnls:
                i = pnl.lpid
                self._adc[i, :] = self.avc[i, :]*pnl.tang
        return self._adc
    @property
    def ada(self):
        if self._ada is None:
            num = len(self.pnls)
            self._ada = zeros((num, 1))
            for pnl in self.pnls:
                i = pnl.lpid
                self._ada[i, 0] = pnl.cdoarea
        return self._ada
    @property
    def cdo(self):
        if self._cdo is None:
            dragarea = 0.0
            for pnl in self.pnls:
                i = pnl.lpid
                dragarea += self.ada[i, 0]
            self._cdo = dragarea/self.sref
        return self._cdo
    @property
    def bvg(self):
        if self._bvg is None:
            num = len(self.strps)
            self._bvg = zeros((num, num))
            for strpi in self.strps:
                i = strpi.lsid
                for strpj in self.strps:
                    j =  strpj.lsid
                    self._bvg[i, j] = strpj.trefftz_velocity(strpi.pnti)*strpi.nrmt
        return self._bvg
    @property
    def bdg(self):
        if self._bdg is None:
            num = len(self.strps)
            self._bdg = zeros((num, num))
            for strp in self.strps:
                i = strp.lsid
                for j in range(num):
                    self._bdg[i, j] = strp.trefftz_drag(self.bvg[i, j])
        return self._bdg
    @property
    def blg(self):
        if self._blg is None:
            num = len(self.strps)
            self._blg = zeros((num, 1))
            for strp in self.strps:
                i = strp.lsid
                self._blg[i, 0] = strp.trefftz_lift()
        return self._blg
    @property
    def byg(self):
        if self._byg is None:
            num = len(self.strps)
            self._byg = zeros((num, 1))
            for strp in self.strps:
                i = strp.lsid
                self._byg[i, 0] = strp.trefftz_yfrc()
        return self._byg
    @property
    def bmg(self):
        if self._bmg is None:
            num = len(self.strps)
            self._bmg = zeros((num, 1))
            for strp in self.strps:
                i = strp.lsid
                self._bmg[i, 0] = strp.pnti.y*self.blg[i, 0]-strp.pnti.z*self.byg[i, 0]
        return self._bmg
    @property
    def bda(self):
        if self._bda is None:
            num = len(self.strps)
            self._bda = zeros((num, 1))
            for strp in self.strps:
                i = strp.lsid
                self._bda[i, 0] = strp.cdoarea
        return self._bda
    @property
    def cdo_ff(self):
        if self._cdo_ff is None:
            dragarea = 0.0
            for strp in self.strps:
                i = strp.lsid
                dragarea += self.bda[i, 0]
            self._cdo_ff = dragarea/self.sref
        return self._cdo_ff
    @property
    def ar(self):
        if self._ar is None:
            self._ar = self.bref**2/self.sref
        return self._ar
    @property
    def lstrpi(self):
        return [strp.lsid for strp in self.strps if strp.msid is None]
    @property
    def mstrpi(self):
        return [strp.lsid for strp in self.strps if strp.msid is not None]
    def set_strip_alpha(self, alpha: list):
        for i, strp in enumerate(self.strps):
            strp._ang = alpha[i]
        self._aic = None
        self._afs = None
        self._gam = None
    @property
    def strip_geometry(self):
        from py2md.classes import MDTable
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('Xle', '.5f')
        table.add_column('Yle', '.5f')
        table.add_column('Zle', '.5f')
        table.add_column('Chord', '.4f')
        table.add_column('Width', '.5f')
        table.add_column('Area', '.6f')
        table.add_column('Dihed', '.4f')
        table.add_column('Incid', '.4f')
        for strp in self.strps:
            j = strp.lsid
            xle = strp.pnti.x
            yle = strp.pnti.y
            zle = strp.pnti.z
            chord = strp.chord
            width = strp.dst
            area = strp.area
            dihed = strp.dihedral
            angle = strp.angle
            table.add_row([j, xle, yle, zle, chord, width, area, dihed, angle])
        return table
    @property
    def panel_geometry(self):
        from py2md.classes import MDTable
        table = MDTable()
        table.add_column('#', 'd')
        table.add_column('X', '.5f')
        table.add_column('Y', '.5f')
        table.add_column('Z', '.5f')
        table.add_column('DX', '.5f')
        for pnl in self.pnls:
            j = pnl.lpid
            x = pnl.pntg.x
            y = pnl.pntg.y
            z = pnl.pntg.z
            dx = pnl.crd
            table.add_row([j, x, y, z, dx])
        return table
    def copy_from_source(self):
        lsys = latticesystem_from_json(self.source, build=False)
        for attr in self.__dict__:
            if attr[0] == '_':
                lsys.__dict__[attr] = self.__dict__[attr].copy()
        lsys.build()
        return lsys
    def __repr__(self):
        return '<LatticeSystem: {:s}>'.format(self.name)
    def __str__(self):
        from py2md.classes import MDTable
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
        table.add_column('# Strips', 'd', data=[len(self.strps)])
        table.add_column('# Panels', 'd', data=[len(self.pnls)])
        table.add_column('# Controls', 'd', data=[len(self.ctrls)])
        outstr += table._repr_markdown_()
        return outstr
    def _repr_markdown_(self):
        return self.__str__()

def latticesystem_from_json(jsonfilepath: str, build: bool=True):
    from .latticesurface import latticesurface_from_json
    from .latticeresult import latticeresult_from_json
    from .latticetrim import latticetrim_from_json
    from pyvlm.tools import masses_from_json, masses_from_data
    from json import load

    with open(jsonfilepath, 'rt') as jsonfile:
        data = load(jsonfile)
    
    from os.path import dirname, join, exists

    path = dirname(jsonfilepath)

    for surfdata in data['surfaces']:
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
    
    name = data['name']
    sfcs = []
    for surfdata in data['surfaces']:
        sfc = latticesurface_from_json(surfdata)
        sfcs.append(sfc)
    bref = data['bref']
    cref = data['cref']
    sref = data['sref']
    xref = data['xref']
    yref = data['yref']
    zref = data['zref']
    rref = Point(xref, yref, zref)
    lsys = LatticeSystem(name, sfcs, bref, cref, sref, rref)

    masses = {}
    if 'masses' in data:
        if isinstance(data['masses'], list):
            masses = masses_from_data(data['masses'])
        elif isinstance(data['masses'], str):
            if data['masses'][-5:] == '.json':
                massfilename = data['masses']
                massfilepath = join(path, massfilename)
            masses = masses_from_json(massfilepath)
    lsys.masses = masses

    if 'cases' in data and build:
        lsys.build()
        for i in range(len(data['cases'])):
            resdata = data['cases'][i]
            if 'trim' in resdata:
                latticetrim_from_json(lsys, resdata)
            else:
                latticeresult_from_json(lsys, resdata)
    
    lsys.source = jsonfilepath

    return lsys

def velocity_matrix(ra: MatrixVector, rb: MatrixVector, rc: MatrixVector):
    
    if ra.shape != rb.shape:
        return ValueError()
    
    numi = rc.shape[0]
    numj = ra.shape[1]
    
    ra = ra.repeat(numi, axis=0)
    rb = rb.repeat(numi, axis=0)
    rc = rc.repeat(numj, axis=1)
    
    a = rc-ra
    b = rc-rb

    am = a.return_magnitude()
    bm = b.return_magnitude()

    # Velocity from Bound Vortex
    adb = elementwise_dot_product(a, b)
    abm = multiply(am, bm)
    dm = multiply(abm, abm+adb)
    axb = elementwise_cross_product(a, b)
    axbm = axb.return_magnitude()
    chki = (axbm == 0.0)
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
