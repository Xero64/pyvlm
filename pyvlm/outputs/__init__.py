from .msh import latticeresult_to_msh
from .md import latticesystem_to_md

def outputs_from_json(sysdct: dict):
    outputs = {}
    for casedct in sysdct['cases']:
        name = casedct['name']
        outputs[name] = []
        if 'outputs' in casedct:
            outputs[name] = outputs[name] + casedct['outputs']
    return outputs
