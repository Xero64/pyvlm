from pyvlm.classes import LatticeSystem


def latticesystem_to_md(lsys: LatticeSystem, mdfilepath: str, outputs: dict=None):
    with open(mdfilepath, 'wt') as mdfile:
        mdfile.write(lsys.__str__())
        for case in lsys.results:
            lres = lsys.results[case]
            mdfile.write('\n')
            mdfile.write(lres.__str__())
            if outputs is not None:
                for output in outputs[case]:
                    if output == 'stability derivatives':
                        mdfile.write('\n')
                        mdfile.write(str(lres.stability_derivatives))
                    if output == 'stability derivatives body':
                        mdfile.write('\n')
                        mdfile.write(str(lres.stability_derivatives_body))
