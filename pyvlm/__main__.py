from json import load
from sys import argv

from .classes import LatticeSystem
from .outputs import latticesystem_to_md, outputs_from_json


def main(jsonfilepath: str = '', mdfilepath: str = ''):

    if jsonfilepath == '':
        if len(argv) == 1:
            print('Specify a .json input file to run and create a .md output file.')
            return
        jsonfilepath = argv[1]

    with open(jsonfilepath, 'rt') as jsonfile:
        sysdct = load(jsonfile)

    sysdct['source'] = jsonfilepath

    sys = LatticeSystem.from_dict(sysdct)

    if mdfilepath == '':
        if len(argv) == 3:
            mdfilepath = argv[2]
        else:
            mdfilepath = jsonfilepath.replace('.json', '.md')

    outputs = outputs_from_json(sysdct)

    latticesystem_to_md(sys, mdfilepath, outputs)
