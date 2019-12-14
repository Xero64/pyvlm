from pyvlm.classes import LatticeSystem

def latticesystem_to_md(lsys: LatticeSystem, mdfilepath: str):
    with open(mdfilepath, 'wt') as mdfile:
        mdfile.write(lsys.__str__())
        for case in lsys.results:
            lres = lsys.results[case]
            mdfile.write('\n')
            mdfile.write(lres.__str__())
