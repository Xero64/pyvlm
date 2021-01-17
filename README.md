# pyvlm
A Vortex Lattice Method implemented in Python for use from terminal or from within Python.

**Installation:**

```
pip install pyvlm
```

**Terminal Usage:**

```
pyvlm aircraft.json
```

The above command generates a "aircraft.md" file containing results.

**JSON Definition File:**

```json
{
    "name": "Test Aircraft",
    "mach": 0.0,
    "sref": 1.358,
    "cref": 0.31,
    "bref": 4.5,
    "xref": 1.175,
    "yref": 0.0,
    "zref": 0.0,
    "surfaces": [
        {
            "name": "Wing",
            "angle": 0.0,
            "mirror": true,
            "numc": 8,
            "cspace": "cosine",
            "xpos": 0.875,
            "ypos": 0.0,
            "zpos": 0.15,
            "sections": [
                {
                    "xpos": 0.0,
                    "ypos": 0.0,
                    "zpos": 0.0,
                    "chord": 0.35,
                    "angle": 0.0,
                    "numb": 3,
                    "bspace": "equal",
                    "airfoil": "NACA 2412",
                    "noload": true
                },
                {
                    "xpos": null,
                    "ypos": 0.21,
                    "zpos": null,
                    "chord": 0.35,
                    "angle": 0.0,
                    "numb": 5,
                    "bspace": "semi-cosine",
                    "airfoil": "NACA 2412"
                },
                {
                    "xpos": null,
                    "ypos": 0.8,
                    "zpos": null,
                    "chord": null,
                    "angle": null,
                    "numb": 20,
                    "bspace": "semi-cosine",
                    "airfoil": "NACA 2412"
                },
                {
                    "xpos": null,
                    "ypos": 1.8,
                    "zpos": null,
                    "chord": null,
                    "angle": null,
                    "numb": 15,
                    "bspace": "cosine",
                    "airfoil": "NACA 2412",
                    "controls": {
                        "aileron": {
                            "xhinge": 0.8,
                            "posgain": 1.0,
                            "neggain": 0.8,
                            "reverse": true,
                            "hvec": {"x": 0.0, "y": 0.0, "z": 0.0}
                        }
                    }
                },
                {
                    "xpos": 0.3,
                    "ypos": 2.25,
                    "zpos": 0.0,
                    "chord": 0.2,
                    "angle": -1.5
                }
            ]
        },
        {
            "name": "Horizontal Tail",
            "mirror": true,
            "numc": 8,
            "cspace": "cosine",
            "xpos": 2.3,
            "ypos": 0.0,
            "zpos": 0.1,
            "angle": -1.8,
            "sections": [
                {
                    "xpos": 0.0,
                    "ypos": 0.0,
                    "zpos": 0.0,
                    "chord": 0.22,
                    "numb": 25,
                    "bspace": "cosine",
                    "controls": {
                        "elevator": {
                            "xhinge": 0.7,
                            "posgain": 1.0,
                            "neggain": 1.0,
                            "reverse": false,
                            "hvec": {"x": 0.0, "y": 0.0, "z": 0.0}
                        }
                    }
                },
                {
                    "xpos": 0.1,
                    "ypos": 0.7,
                    "zpos": 0.0,
                    "chord": 0.18
                }
            ]
        },
        {
            "name": "Vertical Tail",
            "angle": 0.0,
            "numc": 5,
            "cspace": "cosine",
            "xpos": 2.3,
            "ypos": 0.0,
            "zpos": 0.1,
            "sections": [
                {
                    "xpos": 0.0,
                    "ypos": 0.0,
                    "zpos": 0.0,
                    "chord": 0.25,
                    "angle": 0.0,
                    "numb": 15,
                    "bspace": "cosine"
                },
                {
                    "xpos": 0.1,
                    "ypos": 0.0,
                    "zpos": 0.4,
                    "chord": 0.15,
                    "angle": 0.0
                }
            ]
        }
    ],
    "cases": [
        {
            "name": "Positive 1g Cruise",
            "trim": "Looping Trim",
            "density": 0.945,
            "speed": 25.0,
            "mass": 20.0,
            "load factor": 1.0
        },
        {
            "name": "Positive 5g Dive",
            "trim": "Looping Trim",
            "density": 0.945,
            "speed": 50.0,
            "mass": 20.0,
            "load factor": 5.0
        },
        {
            "name": "Negative 3g Dive",
            "trim": "Looping Trim",
            "density": 0.945,
            "speed": 50.0,
            "mass": 20.0,
            "load factor": -3.0
        },
        {
            "name": "60deg Banked Turn Cruise",
            "trim": "Turning Trim",
            "density": 0.945,
            "speed": 25.0,
            "mass": 20.0,
            "bank angle": 60.0
        },
        {
            "name": "Positive 1g Cruise + 15deg Side Slip",
            "inherit": "Positive 1g Cruise",
            "beta": 15.0
        },
        {
            "name": "Positive 1g Cruise + 15deg Elevator",
            "inherit": "Positive 1g Cruise",
            "elevator": 15.0
        },
        {
            "name": "Positive 1g Cruise - 15deg Elevator",
            "inherit": "Positive 1g Cruise",
            "elevator": -15.0
        }
    ]
}
```

**Typical Python Script File "aircraft.py":**

```python
#%% Import Dependencies
from IPython.display import display_markdown
from pyvlm import latticesystem_from_json
from pyvlm.outputs.msh import latticeresult_to_msh
from pyvlm.outputs.prf import latticeresult_to_prf

#%% Import Geometry
jsonfilepath = '../files/Aircraft.json'
lsys = latticesystem_from_json(jsonfilepath)

#%% Display System
display_markdown(lsys)

#%% Display Results
for case in lsys.results:
    lres = lsys.results[case]
    display_markdown(lres)

#%% Mesh File Output
lres = lsys.results['Positive 1g Cruise + 15deg Side Slip']
latticeresult_to_msh(lres, '../results/Aircraft.msh')

#%% Pessure File Output
latticeresult_to_prf(lsys, '../results/Aircraft_pressures.json')

#%% 5g Trim Case
ltrm = lsys.results['Positive 5g Dive']

#%% Plot Lift Distribution
axl = ltrm.plot_trefftz_lift_force_distribution()

#%% Plot Y Force Distribution
axy = ltrm.plot_trefftz_side_force_distribution()

#%% Print Strip Forces
display_markdown(ltrm.strip_forces)

#%% Print Strip Coefficients
display_markdown(ltrm.strip_coefficients)

#%% Print Panel Forces
display_markdown(ltrm.panel_forces)

#%% Print Total Loads
display_markdown(ltrm.surface_loads)
```

**Mesh File Output:**

You can generate a Gmsh mesh file (*.msh) directly from a Python script using the following code snippet.

```python
lres = lsys.results['Positive 1g Cruise + 15deg Side Slip']
latticeresult_to_msh(lres, '../results/Aircraft.msh')
```

This will output a mesh file to the specified location, which can then be viewed in Gmsh. The latest version of Gmsh can be downloaded at:

http://gmsh.info/

Use File > Open in Gmsh to open the mesh file with the pressure results.

A sample of the aircraft shown in Gmsh is captured below. Consult Gmsh help to operate Gmsh.

![](https://github.com/Xero64/pyvlm/raw/main/Readme.png)
