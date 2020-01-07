from json import dump
from pyvlm.classes import LatticeSystem
from numpy.matlib import zeros

def latticeresult_to_prf(lsys: LatticeSystem, prffilepath: str):
    data = {}
    data['name'] = lsys.name
    srfclist = []
    for srfc in lsys.srfcs:
        srfcdict = {}
        xle, yle, zle, chords = [], [], [], []
        csps = srfc.cspace
        cspace = []
        for csp in csps:
            cspace.append(csp[0])
        cspace.append(csp[-1])
        for i in range(srfc.pnts.shape[0]):
            pntle = srfc.pnts[i, 0]
            pntte = srfc.pnts[i, -1]
            x, y, z = pntle.x, pntle.y, pntle.z
            chord = pntte.x-pntle.x
            xle.append(x)
            yle.append(y)
            zle.append(z)
            chords.append(chord)
        srfcdict['name'] = srfc.name
        srfcdict['xle'] = xle
        srfcdict['yle'] = yle
        srfcdict['zle'] = zle
        srfcdict['chords'] = chords
        srfcdict['cspace'] = cspace
        srfclist.append(srfcdict)
    data['surfaces'] = srfclist
    resultlist = []
    for case in lsys.results:
        result = lsys.results[case]
        dynprs = result.qfs
        casedata = {}
        casedata['name'] = case
        casedata['surfaces'] = []
        for srfc in lsys.srfcs:
            prtop = zeros(srfc.pnls.shape)
            prbot = zeros(srfc.pnls.shape)
            for i in range(srfc.pnls.shape[0]):
                for j in range(srfc.pnls.shape[1]):
                    pnl = srfc.pnls[i, j]
                    lpid = pnl.lpid
                    pnlfrc = result.nfres.nffrc[lpid, 0]
                    frcnrm = pnl.nrml*pnlfrc
                    deltap = frcnrm/pnl.area
                    prtop[i, j] = dynprs+deltap/2
                    prbot[i, j] = dynprs-deltap/2
            srfresdata = {}
            srfresdata['name'] = srfc.name
            srfresdata['top pressure [Pa]'] = prtop.tolist()
            srfresdata['bottom pressure [Pa]'] = prbot.tolist()
            casedata['surfaces'].append(srfresdata)
        resultlist.append(casedata)
    data['results'] = resultlist
    with open(prffilepath, 'wt') as prffile:
        dump(data, prffile, indent=4)
