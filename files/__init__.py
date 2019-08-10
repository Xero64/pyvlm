from inspect import getfile, currentframe
from os.path import abspath, join, basename
from glob import glob

path = abspath(getfile(currentframe())).replace('__init__.py', '')

query = join(path, '*.json')
filepaths = glob(query)

filedict = {basename(filepath): filepath for filepath in filepaths}

def load_package_file(filename):
    from ..classes.latticesystem import latticesystem_from_json
    filepath = filedict[filename]
    return latticesystem_from_json(filepath)
