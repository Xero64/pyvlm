from sys import argv

def main():
    if len(argv) == 1:
        print('Specify a .json input file to run and create a .md output file.')
        return
    
    from pyvlm import latticesystem_from_json
    from pyvlm.outputs import latticesystem_to_md
    
    jsonfilepath = argv[1]

    lsys = latticesystem_from_json(jsonfilepath)

    if len(argv) == 3:
        mdfilepath = argv[2]
    else:
        mdfilepath = jsonfilepath.replace('.json', '.md')

    latticesystem_to_md(lsys, mdfilepath)
