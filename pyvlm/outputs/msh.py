from pyvlm.classes import LatticeResult, LatticeSurface

class SurfaceMesh(object):
    sid = None
    srfc = None
    pnts = None
    pnls = None
    def __init__(self, srfc: LatticeSurface):
        self.srfc = srfc
        self.pnts = []
        for i in range(srfc.pnts.shape[0]):
            for j in range(srfc.pnts.shape[1]):
                self.pnts.append(srfc.pnts[i, j])
        self.pnls = []
        for i in range(srfc.pnls.shape[0]):
            for j in range(srfc.pnls.shape[1]):
                self.pnls.append(srfc.pnls[i, j])
    def number_points(self, pntid: int):
        for pnt in self.pnts:
            pnt.pntid = pntid
            pntid += 1
        return pntid
    def min_max_coordinates(self):
        x = [pnt.x for pnt in self.pnts]
        y = [pnt.y for pnt in self.pnts]
        z = [pnt.z for pnt in self.pnts]
        return min(x), min(y), min(z), max(x), max(y), max(z)

def latticeresult_to_msh(lres: LatticeResult, mshfilepath: str):
    pnts = []
    pntid = 1
    sid = 1
    srfcmshs = [SurfaceMesh(srfc) for srfc in lres.sys.srfcs]
    for srfcmsh in srfcmshs:
        srfcmsh.sid = sid
        sid += 1
        pnts += srfcmsh.pnts
        pntid = srfcmsh.number_points(pntid)
    ns = len(srfcmshs)
    pnlids = []
    for pnl in lres.sys.pnls:
        pnlids.append(pnl.lpid+1)
    minpnt = 1
    maxpnt = pntid-1
    lenpnt = pntid-1
    minpnl = min(pnlids)
    maxpnl = max(pnlids)
    lenpnl = len(pnlids)
    with open(mshfilepath, 'wt') as mshfile:
        mshfile.write('$MeshFormat\n')
        mshfile.write('4.1 0 8\n')
        mshfile.write('$EndMeshFormat\n')
        mshfile.write('$PartitionedEntities\n')
        mshfile.write('{:d}\n'.format(ns))
        mshfile.write('0\n')
        mshfile.write('0 0 {:d} 0\n'.format(ns))
        for srfcmsh in srfcmshs:
            mshfile.write('{:d} 2 1 1 {:d} '.format(srfcmsh.sid, srfcmsh.sid))
            mshfile.write('{:} {:} {:} {:} {:} {:} 0 0\n'.format(*srfcmsh.min_max_coordinates()))
        mshfile.write('$EndPartitionedEntities\n')
        mshfile.write('$Nodes\n')
        mshfile.write('{:d} {:d} {:d} {:d}\n'.format(ns, lenpnt, minpnt, maxpnt))
        for srfcmsh in srfcmshs:
            sid = srfcmsh.sid
            nsp = len(srfcmsh.pnts)
            mshfile.write('{:d} {:d} {:d} {:d}\n'.format(2, sid, 0, nsp))
            frmstr = '{:d}\n'
            for pnt in srfcmsh.pnts:
                mshfile.write(frmstr.format(pnt.pntid))
            frmstr = '{:} {:} {:}\n'
            for pnt in srfcmsh.pnts:
                x, y, z = pnt.x, pnt.y, pnt.z
                mshfile.write(frmstr.format(x, y, z))
        mshfile.write('$EndNodes\n')
        mshfile.write('$Elements\n')
        mshfile.write('{:d} {:d} {:d} {:d}\n'.format(ns, lenpnl, minpnl, maxpnl))
        for srfcmsh in srfcmshs:
            sid = srfcmsh.sid
            nsp = len(srfcmsh.pnls)
            mshfile.write('{:d} {:d} {:d} {:d}\n'.format(2, sid, 3, nsp))
            for pnl in srfcmsh.pnls:
                outstr = '{:d}'.format(pnl.lpid+1)
                outstr += ' {:d}'.format(pnl.pnts[0].pntid)
                outstr += ' {:d}'.format(pnl.pnts[1].pntid)
                outstr += ' {:d}'.format(pnl.pnts[3].pntid)
                outstr += ' {:d}'.format(pnl.pnts[2].pntid)
                outstr += '\n'
                mshfile.write(outstr)
        mshfile.write('$EndElements\n')
        optstr = ''
        optstr += 'Mesh.Lines = 0;\n'
        optstr += 'Mesh.LineNumbers = 0;\n'
        optstr += 'Mesh.SurfaceEdges = 1;\n'
        optstr += 'Mesh.SurfaceFaces = 0;\n'
        optstr += 'Mesh.SurfaceNumbers = 0;\n'
        optstr += 'Mesh.VolumeEdges = 0;\n'
        optstr += 'Mesh.VolumeFaces = 0;\n'
        optstr += 'Mesh.VolumeNumbers = 0;\n'
        view = 0
        mshfile.write('$ElementData\n')
        mshfile.write('1\n')
        mshfile.write('"Gamma [m/s]"\n')
        mshfile.write('1\n')
        mshfile.write('0.0\n')
        mshfile.write('3\n')
        mshfile.write('0\n')
        mshfile.write('1\n')
        mshfile.write('{:d}\n'.format(lenpnl))
        frmstr = '{:d} {:f}\n'
        for pnl in lres.sys.pnls:
            gamma = lres.nfres.gamma[pnl.lpid, 0]
            mshfile.write(frmstr.format(pnl.lpid+1, gamma))
        mshfile.write('$EndElementData\n')
        optstr += 'View[{:d}].Light = 0;\n'.format(view)
        optstr += 'View[{:d}].SaturateValues = 1;\n'.format(view)
        optstr += 'View[{:d}].Visible = 0;\n'.format(view)
        view += 1
        mshfile.write('$ElementData\n')
        mshfile.write('1\n')
        mshfile.write('"Pressure Delta [Pa]"\n')
        mshfile.write('1\n')
        mshfile.write('0.0\n')
        mshfile.write('3\n')
        mshfile.write('0\n')
        mshfile.write('1\n')
        mshfile.write('{:d}\n'.format(lenpnl))
        frmstr = '{:d} {:f}\n'
        for pnl in lres.sys.pnls:
            dp = lres.nfres.nffrc[pnl.lpid, 0]*pnl.nrml/pnl.area
            mshfile.write(frmstr.format(pnl.lpid+1, dp))
        mshfile.write('$EndElementData\n')
        optstr += 'View[{:d}].Light = 0;\n'.format(view)
        optstr += 'View[{:d}].SaturateValues = 1;\n'.format(view)
        optstr += 'View[{:d}].Visible = 0;\n'.format(view)
        view += 1
        p = lres.qfs
        mshfile.write('$ElementData\n')
        mshfile.write('1\n')
        mshfile.write('"Pressure Bottom [Pa]"\n')
        mshfile.write('1\n')
        mshfile.write('0.0\n')
        mshfile.write('3\n')
        mshfile.write('0\n')
        mshfile.write('1\n')
        mshfile.write('{:d}\n'.format(lenpnl))
        frmstr = '{:d} {:f}\n'
        for pnl in lres.sys.pnls:
            dp = lres.nfres.nffrc[pnl.lpid, 0]*pnl.nrml/pnl.area
            mshfile.write(frmstr.format(pnl.lpid+1, p-dp/2))
        mshfile.write('$EndElementData\n')
        optstr += 'View[{:d}].Light = 0;\n'.format(view)
        optstr += 'View[{:d}].SaturateValues = 1;\n'.format(view)
        optstr += 'View[{:d}].Visible = 0;\n'.format(view)
        view += 1
        p = lres.qfs
        mshfile.write('$ElementData\n')
        mshfile.write('1\n')
        mshfile.write('"Pressure Top [Pa]"\n')
        mshfile.write('1\n')
        mshfile.write('0.0\n')
        mshfile.write('3\n')
        mshfile.write('0\n')
        mshfile.write('1\n')
        mshfile.write('{:d}\n'.format(lenpnl))
        frmstr = '{:d} {:f}\n'
        for pnl in lres.sys.pnls:
            dp = lres.nfres.nffrc[pnl.lpid, 0]*pnl.nrml/pnl.area
            mshfile.write(frmstr.format(pnl.lpid+1, p+dp/2))
        mshfile.write('$EndElementData\n')
        optstr += 'View[{:d}].Light = 0;\n'.format(view)
        optstr += 'View[{:d}].SaturateValues = 1;\n'.format(view)
        optstr += 'View[{:d}].Visible = 0;\n'.format(view)
        view += 1
    optfilepath = mshfilepath + '.opt'
    with open(optfilepath, 'wt') as optfile:
        optfile.write(optstr)
