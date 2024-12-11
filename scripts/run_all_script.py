#%%
# Import Dependencies
from os import curdir
from os.path import join
from subprocess import run
from pyvlm.__main__ import main

#%%
# Run All Test Scripts
test_list = [
    'Test_straight',
    'Test_taper',
    'Test_sweep',
    'Test_dihedral',
    'Test_twist',
    'Test_camber',
    'Test_twist_taper'
]

for test in test_list:
    jsonfilepath = join(curdir, '..', 'files', test + '.json')
    mdfilepath = join(curdir, '..', 'results', test + '.md')
    print(f'Running {test}...')
    main(jsonfilepath, mdfilepath)
