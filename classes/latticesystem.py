from time import time
from math import pi, atan2
from numpy.matlib import zeros, empty
from numpy.linalg import solve, inv
from pygeom.geom3d import Point, Vector, ihat, jhat, zero_vector

class LatticeSystem(object):
    name = None
    srfcs = None
    strps = None
    mstrp = None
    pnls = None
    mpnl = None
    bref = None
    cref = None
    sref = None
    rref = None
    ctrls = None
    _gam = None
    _avg = None
    _afg = None
    _amg = None
    _avc = None
    _aic = None
    _afs = None
    _adc = None
    _bvg = None
    _bdg = None
    _blg = None
    _byg = None
    _bmg = None
    _ar = None
    _strpy = None
    def __init__(self, name: str, srfcs: list, bref: float, cref: float, sref: float, rref: Point):
        self.name = name
        self.srfcs = srfcs
        self.bref = bref
        self.cref = cref
        self.sref = sref
        self.rref = rref
        self.mesh()
        self.inherit()
        self.build()
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
        ind = 2
        for srfc in self.srfcs:
            for sht in srfc.shts:
                for control in sht.ctrls:
                    if control not in self.ctrls:
                        # self.ctrls[control] = (ind, ind+1, ind+2, ind+3, ind+4, ind+5)
                        self.ctrls[control] = (ind, ind+1, ind+2, ind+3)
                        ind += 2
    def build(self):
        self.avc
        self.aic
        self.afs
        self.avg
        self.afg
        self.amg
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
    def gam(self):
        if self._gam is None:
            num = len(self.pnls)
            numc = len(self.ctrls)
            gamma = -solve(self.aic, self.afs)
            self._gam = empty((num, 2+4*numc), dtype=Vector)
            for i in range(num):
                for j in range(2+4*numc):
                    self._gam[i, j] = Vector(gamma[i, j*3], gamma[i, j*3+1], gamma[i, j*3+2])
            # self._gam = -solve(self.aic, self.afs)
        return self._gam
    @property
    def avg(self):
        if self._avg is None:
            print('Building Induced Velocity Matrix.')
            start = time()
            num = len(self.pnls)
            self._avg = empty((num, num), dtype=Vector)
            for pnli in self.pnls:
                i = pnli.lpid
                for pnlj in self.pnls:
                    j = pnlj.lpid
                    self._avg[i, j] = pnlj.velocity(pnli.pnti)
            finish = time()
            elapsed = finish-start
            print(f'Built Induced Velocity Matrix in {elapsed:.3f} seconds.')
        return self._avg
    @property
    def afg(self):
        if self._afg is None:
            print('Building Induced Force Matrix.')
            start = time()
            num = len(self.pnls)
            self._afg = empty((num, num), dtype=Vector)
            for pnl in self.pnls:
                i = pnl.lpid
                for j in range(num):
                    self._afg[i, j] = self.avg[i, j]**pnl.leni
            finish = time()
            elapsed = finish-start
            print(f'Built Induced Force Matrix in {elapsed:.3f} seconds.')
        return self._afg
    @property
    def amg(self):
        if self._amg is None:
            print('Building Induced Moment Matrix.')
            start = time()
            num = len(self.pnls)
            self._amg = empty((num, num), dtype=Vector)
            for pnl in self.pnls:
                i = pnl.lpid
                for j in range(num):
                    self._amg[i, j] = (pnl.pnti-self.rref)**self.afg[i, j]
            finish = time()
            elapsed = finish-start
            print(f'Built Induced Moment Matrix in {elapsed:.3f} seconds.')
        return self._amg
    @property
    def avc(self):
        if self._avc is None:
            print('Building Control Velocity Matrix.')
            start = time()
            num = len(self.pnls)
            self._avc = empty((num, num), dtype=Vector)
            for pnli in self.pnls:
                i = pnli.lpid
                for pnlj in self.pnls:
                    j = pnlj.lpid
                    self._avc[i, j] = pnlj.velocity(pnli.pntc)
            finish = time()
            elapsed = finish-start
            print(f'Built Control Velocity Matrix in {elapsed:.3f} seconds.')
        return self._avc
    @property
    def aic(self):
        if self._aic is None:
            print('Building AIC Matrix.')
            start = time()
            num = len(self.pnls)
            self._aic = zeros((num, num))
            for pnl in self.pnls:
                i = pnl.lpid
                for j in range(num):
                    self._aic[i, j] = self.avc[i, j]*pnl.nrml
            finish = time()
            elapsed = finish-start
            print(f'Built AIC Matrix in {elapsed:.3f} seconds.')
        return self._aic
    @property
    def adc(self):
        if self._adc is None:
            num = len(self.pnls)
            self._adc = zeros((num, num))
            for pnl in self.pnls:
                i = pnl.lpid
                for j in range(num):
                    self._adc[i, j] = self.avc[i, j]*pnl.tang
        return self._adc
    @property
    def afs(self):
        if self._afs is None:
            print('Building Freestream and Control Matrix.')
            start = time()
            num = len(self.pnls)
            numc = len(self.ctrls)
            self._afs = zeros((num, 6+12*numc))
            for pnl in self.pnls:
                lpid = pnl.lpid
                rpx, rpy, rpz = pnl.pntc.x, pnl.pntc.y, pnl.pntc.z
                nmx, nmy, nmz = pnl.nrml.x, pnl.nrml.y, pnl.nrml.z
                self._afs[lpid, 0] = nmx
                self._afs[lpid, 1] = nmy
                self._afs[lpid, 2] = nmz
                self._afs[lpid, 3] = nmz*rpy - nmy*rpz
                self._afs[lpid, 4] = nmx*rpz - nmz*rpx
                self._afs[lpid, 5] = nmy*rpx - nmx*rpy
            for srfc in self.srfcs:
                for sht in srfc.shts:
                    for control in sht.ctrls:
                        ctrl = sht.ctrls[control]
                        ctup = self.ctrls[control]
                        for pnl in ctrl.pnls:
                            lpid = pnl.lpid
                            dndlp = pnl.dndl(ctrl.posgain, ctrl.uhvec)
                            nmx, nmy, nmz = dndlp.x, dndlp.y, dndlp.z
                            self._afs[lpid, ctup[0]*3+0] = nmx
                            self._afs[lpid, ctup[0]*3+1] = nmy
                            self._afs[lpid, ctup[0]*3+2] = nmz
                            self._afs[lpid, ctup[1]*3+0] = nmz*rpy - nmy*rpz
                            self._afs[lpid, ctup[1]*3+1] = nmx*rpz - nmz*rpx
                            self._afs[lpid, ctup[1]*3+2] = nmy*rpx - nmx*rpy
                            dndln = pnl.dndl(ctrl.neggain, ctrl.uhvec)
                            nmx, nmy, nmz = dndln.x, dndln.y, dndln.z
                            self._afs[lpid, ctup[2]*3+0] = nmx
                            self._afs[lpid, ctup[2]*3+1] = nmy
                            self._afs[lpid, ctup[2]*3+2] = nmz
                            self._afs[lpid, ctup[3]*3+0] = nmz*rpy - nmy*rpz
                            self._afs[lpid, ctup[3]*3+1] = nmx*rpz - nmz*rpx
                            self._afs[lpid, ctup[3]*3+2] = nmy*rpx - nmx*rpy
            finish = time()
            elapsed = finish-start
            print(f'Built Freestream and Control Matrix in {elapsed:.3f} seconds.')
        return self._afs
    @property
    def bvg(self):
        if self._bvg is None:
            print('Building Trefftz Induced Velocity Matrix.')
            start = time()
            num = len(self.strps)
            self._bvg = zeros((num, num))
            for strpi in self.strps:
                i = strpi.lsid
                for strpj in self.strps:
                    j =  strpj.lsid
                    self._bvg[i, j] = strpj.trefftz_velocity(strpi.pnti)*strpi.nrmt
            finish = time()
            elapsed = finish-start
            print(f'Built Trefftz Induced Velocity Matrix in {elapsed:.3f} seconds.')
        return self._bvg
    @property
    def bdg(self):
        if self._bdg is None:
            print('Building Trefftz Induced Drag Matrix.')
            start = time()
            num = len(self.strps)
            self._bdg = zeros((num, num))
            for strp in self.strps:
                i = strp.lsid
                for j in range(num):
                    self._bdg[i, j] = -strp.dst*self.bvg[i, j]/2
            finish = time()
            elapsed = finish-start
            print(f'Built Trefftz Induced Drag Matrix in {elapsed:.3f} seconds.')
        return self._bdg
    @property
    def blg(self):
        if self._blg is None:
            print('Building Trefftz Induced Lift Matrix.')
            start = time()
            num = len(self.strps)
            self._blg = zeros((num, 1))
            for strp in self.strps:
                i = strp.lsid
                self._blg[i, 0] = strp.lent.y
            finish = time()
            elapsed = finish-start
            print(f'Built Trefftz Induced Lift Matrix in {elapsed:.3f} seconds.')
        return self._blg
    @property
    def byg(self):
        if self._byg is None:
            print('Building Trefftz Induced Y-Force Matrix.')
            start = time()
            num = len(self.strps)
            self._byg = zeros((num, 1))
            for strp in self.strps:
                i = strp.lsid
                self._byg[i, 0] = -strp.lent.z
            finish = time()
            elapsed = finish-start
            print(f'Built Trefftz Induced Y-Force Matrix in {elapsed:.3f} seconds.')
        return self._byg
    @property
    def bmg(self):
        if self._bmg is None:
            print('Building Trefftz Induced X-Moment Matrix.')
            start = time()
            num = len(self.strps)
            self._bmg = zeros((num, 1))
            for strp in self.strps:
                i = strp.lsid
                self._bmg[i, 0] = strp.pnti.y*self.blg[i, 0]-strp.pnti.z*self.byg[i, 0]
            finish = time()
            elapsed = finish-start
            print(f'Built Trefftz Induced X-Moment Matrix in {elapsed:.3f} seconds.')
        return self._bmg
    @property
    def ar(self):
        if self._ar is None:
            self._ar = self.bref**2/self.sref
        return self._ar
    def set_strip_alpha(self, alpha: list):
        for i, strp in enumerate(self.strps):
            strp._ang = alpha[i]
        self.reset()
    def print_strip_geometry(self, filepath: str=''):
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
        if filepath != '':
            with open(filepath, 'wt') as resfile:
                resfile.write(table._repr_markdown_())
        else:
            print(table._repr_markdown_())
    def print_panel_geometry(self, filepath: str=''):
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
        if filepath != '':
            with open(filepath, 'wt') as resfile:
                resfile.write(table._repr_markdown_())
        else:
            print(table._repr_markdown_())
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

def latticesystem_from_json(jsonfilepath: str):
    from .latticesurface import latticesurface_from_json
    from json import load
    with open(jsonfilepath, 'rt') as jsonfile:
        data = load(jsonfile)
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
    sys = LatticeSystem(name, sfcs, bref, cref, sref, rref)
    return sys
