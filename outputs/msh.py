from pyvlm.classes import LatticeResult

def latticeresult_to_msh(lres: LatticeResult, mshfilepath: str):
    # Add point ids
    pnts = []
    pntid = 1
    pntids = []
    for srfc in lres.sys.srfcs:
        for i in range(srfc.pnts.shape[0]):
            for j in range(srfc.pnts.shape[1]):
                pnts.append(srfc.pnts[i, j]) 
                srfc.pnts[i, j].pntid = pntid
                pntids.append(pntid)
                pntid += 1
    for srfc in lres.sys.srfcs:
        for i in range(srfc.pnls.shape[0]):
            for j in range(srfc.pnls.shape[1]):
                pnl = srfc.pnls[i, j]
                pnts.append(pnl.pntc) 
                pnl.pntc.pntid = pntid
                pntids.append(pntid)
                pntid += 1
    for srfc in lres.sys.srfcs:
        for i in range(srfc.pnls.shape[0]):
            for j in range(srfc.pnls.shape[1]):
                pnl = srfc.pnls[i, j]
                pnts.append(pnl.pnti) 
                pnl.pnti.pntid = pntid
                pntids.append(pntid)
                pntid += 1
    pnlids = []
    for pnl in lres.sys.pnls:
        pnlids.append(pnl.lpid+1)
    minpnt = min(pntids)
    maxpnt = max(pntids)
    lenpnt = len(pntids)
    minpnl = min(pnlids)
    maxpnl = max(pnlids)
    lenpnl = len(pnlids)
    with open(mshfilepath, 'wt') as mshfile:
        mshfile.write('$MeshFormat\n')
        mshfile.write('4.1 0 8\n')
        mshfile.write('$EndMeshFormat\n')
        mshfile.write('$Nodes\n')
        mshfile.write('{:d} {:d} {:d} {:d}\n'.format(1, lenpnt, minpnt, maxpnt))
        mshfile.write('{:d} {:d} {:d} {:d}\n'.format(2, 1, 0, lenpnt))
        frmstr = '{:d}\n'
        for pntid in pntids:
            mshfile.write(frmstr.format(pntid))
        frmstr = '{:} {:} {:}\n'
        for pnt in pnts:
            x = pnt.x
            y = pnt.y
            z = pnt.z
            mshfile.write(frmstr.format(x, y, z))
        mshfile.write('$EndNodes\n')
        mshfile.write('$Elements\n')
        mshfile.write('{:d} {:d} {:d} {:d}\n'.format(1, lenpnl, minpnl, maxpnl))
        mshfile.write('{:d} {:d} {:d} {:d}\n'.format(2, 1, 3, lenpnl))
        for pnl in lres.sys.pnls:
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
        optstr += 'Mesh.SurfaceEdges = 0;\n'
        optstr += 'Mesh.SurfaceFaces = 0;\n'
        optstr += 'Mesh.SurfaceNumbers = 0;\n'
        optstr += 'Mesh.VolumeEdges = 0;\n'
        optstr += 'Mesh.VolumeFaces = 0;\n'
        optstr += 'Mesh.VolumeNumbers = 0;\n'
        view = 0
        mshfile.write('$NodeData\n')
        mshfile.write('1\n')
        mshfile.write('"Panel Normal"\n')
        mshfile.write('1\n')
        mshfile.write('0.0\n')
        mshfile.write('3\n')
        mshfile.write('0\n')
        mshfile.write('3\n')
        mshfile.write('{:d}\n'.format(lenpnl))
        frmstr = '{:d} {:f} {:f} {:f}\n'
        for pnl in lres.sys.pnls:
            pntid = pnl.pntc.pntid
            nrml = pnl.nrml
            mshfile.write(frmstr.format(pntid, nrml.x, nrml.y, nrml.z))
        mshfile.write('$EndNodeData\n')
        view += 1
        mshfile.write('$ElementData\n')
        mshfile.write('1\n')
        mshfile.write('"Gamma"\n')
        mshfile.write('1\n')
        mshfile.write('0.0\n')
        mshfile.write('3\n')
        mshfile.write('0\n')
        mshfile.write('1\n')
        mshfile.write('{:d}\n'.format(lenpnl))
        frmstr = '{:d} {:f}\n'
        for pnl in lres.sys.pnls:
            gamma = lres.gam[pnl.lpid]
            mshfile.write(frmstr.format(pnl.lpid+1, gamma))
        mshfile.write('$EndElementData\n')
        # optstr += 'View[{:d}].CustomMax = 1;\n'.format(view)
        # optstr += 'View[{:d}].CustomMin = 0;\n'.format(view)
        optstr += 'View[{:d}].Light = 0;\n'.format(view)
        # optstr += 'View[{:d}].RangeType = 2;\n'.format(view)
        optstr += 'View[{:d}].SaturateValues = 1;\n'.format(view)
        optstr += 'View[{:d}].Visible = 0;\n'.format(view)
        view += 1
        mshfile.write('$NodeData\n')
        mshfile.write('1\n')
        mshfile.write('"Panel Forces"\n')
        mshfile.write('1\n')
        mshfile.write('0.0\n')
        mshfile.write('3\n')
        mshfile.write('0\n')
        mshfile.write('3\n')
        mshfile.write('{:d}\n'.format(lenpnl))
        frmstr = '{:d} {:f} {:f} {:f}\n'
        for pnl in lres.sys.pnls:
            pntid = pnl.pnti.pntid
            nffrc = lres.nffrc[pnl.lpid]
            mshfile.write(frmstr.format(pntid, nffrc.x, nffrc.y, nffrc.z))
        mshfile.write('$EndNodeData\n')
        view += 1
    optfilepath = mshfilepath + '.opt'
    with open(optfilepath, 'wt') as optfile:
        optfile.write(optstr)
    #     for lcname in group.loadcases:
    #         for crit in group.result:
    #             lcres = group.result[crit][lcname]
    #             mshfile.write('$ElementData\n')
    #             mshfile.write('1\n')
    #             mshfile.write('"{:s} - {:s}"\n'.format(lcname, crit))
    #             mshfile.write('1\n')
    #             mshfile.write('0.0\n')
    #             mshfile.write('3\n')
    #             mshfile.write('0\n')
    #             mshfile.write('1\n')
    #             rescnt = 0
    #             for eid in lcres:
    #                 if eid in eidlst:
    #                     rescnt += 1
    #             mshfile.write('{:d}\n'.format(rescnt))
    #             frmstr = '{:d} {:}\n'
    #             for eid in lcres:
    #                 if eid in eidlst:
    #                     mshfile.write(frmstr.format(eid, 1/lcres[eid].RF))
    #             mshfile.write('$EndElementData\n')
    #             optstr += 'View[{:d}].CustomMax = 1;\n'.format(view)
    #             optstr += 'View[{:d}].CustomMin = 0;\n'.format(view)
    #             optstr += 'View[{:d}].Light = 0;\n'.format(view)
    #             optstr += 'View[{:d}].RangeType = 2;\n'.format(view)
    #             optstr += 'View[{:d}].SaturateValues = 1;\n'.format(view)
    #             optstr += 'View[{:d}].Visible = 0;\n'.format(view)
    #             view += 1
    # optfilepath = mshfilepath + '.opt'
    # with open(optfilepath, 'wt') as optfile:
    #     optfile.write(optstr)

# def msh_float(value: float):
#     if value.is_integer():
#         return '{:}.0'.format(value)
#     else:
#         return '{:}'.format(value)
