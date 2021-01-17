from numpy.matlib import zeros

class Mass(object):
    name = None
    mass = None
    xcm = None
    ycm = None
    zcm = None
    Ixx = None
    Ixy = None
    Ixz = None
    Iyy = None
    Iyz = None
    Izz = None
    def __init__(self, name: str, mass: float=0.0,
                 xcm: float=0.0, ycm: float=0.0, zcm: float=0.0,
                 Ixx: float=0.0, Ixy: float=0.0, Ixz: float=0.0,
                 Iyy: float=0.0, Iyz: float=0.0, Izz: float=0.0):
        self.name = name
        self.mass = mass
        self.xcm = xcm
        self.ycm = ycm
        self.zcm = zcm
        self.Ixx = Ixx
        self.Ixy = Ixy
        self.Ixz = Ixz
        self.Iyy = Iyy
        self.Iyz = Iyz
        self.Izz = Izz
    def mass_matrix(self):
        M = zeros((3, 3))
        M[0, 0] = self.mass
        M[1, 1] = self.mass
        M[2, 2] = self.mass
        return M
    def moment_of_inertia_matrix(self):
        I = zeros((3, 3))
        I[0, 0] = self.Ixx
        I[1, 1] = self.Iyy
        I[2, 2] = self.Izz
        I[0, 1] = self.Ixy
        I[1, 0] = self.Ixy
        I[0, 2] = self.Ixz
        I[2, 0] = self.Ixz
        I[1, 2] = self.Iyz
        I[2, 1] = self.Iyz
        return I
    def __repr__(self):
        return '<Mass {:s}>'.format(self.name)
    def __str__(self):
        frmstr = 'Mass\t{:s}\t{:}'
        return frmstr.format(self.name, self.mass)
    def __format__(self, format_spec: str):
        frmstr = 'Mass\t{:s}\t{:'
        frmstr += format_spec
        frmstr += '}'
        return frmstr.format(self.name, self.mass)

class MassCollection(Mass):
    masses = None
    def __init__(self, name: str, masses: dict):
        super(MassCollection, self).__init__(name)
        self.masses = masses
        self.update()
    def update(self):
        self.mass = 0.0
        mx = 0.0
        my = 0.0
        mz = 0.0
        for m in self.masses:
            mass = self.masses[m]
            self.mass += mass.mass
            mx += mass.mass*mass.xcm
            my += mass.mass*mass.ycm
            mz += mass.mass*mass.zcm
        self.xcm = mx/self.mass
        self.ycm = my/self.mass
        self.zcm = mz/self.mass
        self.Ixx = 0.0
        self.Ixy = 0.0
        self.Ixz = 0.0
        self.Iyy = 0.0
        self.Iyz = 0.0
        self.Izz = 0.0
        for m in self.masses:
            mass = self.masses[m]
            x = mass.xcm-self.xcm
            y = mass.ycm-self.ycm
            z = mass.zcm-self.zcm
            self.Ixx += mass.Ixx + mass.mass*(y**2 + z**2)
            self.Iyy += mass.Iyy + mass.mass*(z**2 + x**2)
            self.Izz += mass.Izz + mass.mass*(x**2 + y**2)
            self.Ixy += mass.Ixy - mass.mass*x*y
            self.Ixz += mass.Ixz - mass.mass*x*z
            self.Iyz += mass.Iyz - mass.mass*y*z
    def __repr__(self):
        return '<MassCollection {:s}>'.format(self.name)
    def __str__(self):
        frmstr = 'MassCollection\t{:s}\t{:}'
        return frmstr.format(self.name, self.mass)
    def __format__(self, format_spec: str):
        frmstr = 'MassCollection\t{:s}\t{:'
        frmstr += format_spec
        frmstr += '}'
        return frmstr.format(self.name, self.mass)

def mass_from_data(massdata: dict):
    name = massdata['name']
    m = massdata['mass']
    xcm = massdata['xcm']
    ycm = massdata['ycm']
    zcm = massdata['zcm']
    mass = Mass(name, m, xcm, ycm, zcm)
    if 'Ixx' in massdata:
        mass.Ixx = massdata['Ixx']
    if 'Ixy' in massdata:
        mass.Ixx = massdata['Ixy']
    if 'Ixz' in massdata:
        mass.Ixx = massdata['Ixz']
    if 'Iyy' in massdata:
        mass.Ixx = massdata['Iyy']
    if 'Iyz' in massdata:
        mass.Ixx = massdata['Iyz']
    if 'Izz' in massdata:
        mass.Ixx = massdata['Izz']
    return mass

def masscol_from_data(masscoldata: dict, masses: dict):
    name = masscoldata['name']
    masslist = masscoldata['masses']
    massdict = {}
    for m in masslist:
        if m in masses:
            massdict[m] = masses[m]
        else:
            print(f'Could not find {m:s} in masses.')
    return MassCollection(name, massdict)

def masses_from_data(massesdata: list):
    masses = {}
    for mdata in massesdata:
        if 'masses' in mdata:
            masscol = masscol_from_data(mdata, masses)
            masses[masscol.name] = masscol
        else:
            mass = mass_from_data(mdata)
            masses[mass.name] = mass
    return masses

def masses_from_json(jsonfilepath: str):
    from json import load

    with open(jsonfilepath, 'rt') as jsonfile:
        massesdata = load(jsonfile)

    return masses_from_data(massesdata)

def mass_table(masses):
    from py2md.classes import MDTable
    masslst = []
    if isinstance(masses, dict):
        masslst = [masses[m] for m in masses]
    elif isinstance(masses, list):
        masslst = masses
    table = MDTable()
    table.add_column('Name', 's')
    table.add_column('Mass', '.3f')
    table.add_column('xcm', '.3f')
    table.add_column('ycm', '.3f')
    table.add_column('zcm', '.3f')
    table.add_column('Ixx', '.1f')
    table.add_column('Ixy', '.1f')
    table.add_column('Ixz', '.1f')
    table.add_column('Iyy', '.1f')
    table.add_column('Iyz', '.1f')
    table.add_column('Izz', '.1f')
    for mass in masslst:
        table.add_row([mass.name, mass.mass,
                       mass.xcm, mass.ycm, mass.zcm,
                       mass.Ixx, mass.Ixy, mass.Ixz,
                       mass.Iyy, mass.Iyz, mass.Izz])
    return table
