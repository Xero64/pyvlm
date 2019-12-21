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

The above command generates a aircraft.md containing results.

**JSON Definition File:**

```json
{
    "name": "Test Aircraft",
    "mach": 0.0,
    "sref": 1.358,
    "cref": 0.31,
    "bref": 4.5,
    "xref": 0.96,
    "yref": 0.0,
    "zref": 0.0,
    "surfaces": [
        {
            "name": "Wing",
            "angle": 0.0,
            "mirror": true,
            "numc": 8,
            "cspace": "cosine",
            "sections": [
                {
                    "xle": 0.875,
                    "yle": 0.0,
                    "zle": 0.15,
                    "chord": 0.35,
                    "angle": 0.0,
                    "numb": 3,
                    "bspace": "equal",
                    "airfoil": "NACA 2412",
                    "noload": true
                },
                {
                    "xle": 0.875,
                    "yle": 0.21,
                    "zle": 0.15,
                    "chord": 0.35,
                    "angle": 0.0,
                    "numb": 5,
                    "bspace": "semi-cosine",
                    "airfoil": "NACA 2412"
                },
                {
                    "xle": 0.875,
                    "yle": 0.8,
                    "zle": 0.15,
                    "chord": 0.35,
                    "angle": 0.0,
                    "numb": 20,
                    "bspace": "semi-cosine",
                    "airfoil": "NACA 2412"
                },
                {
                    "xle": 0.875,
                    "yle": 1.8,
                    "zle": 0.15,
                    "chord": 0.25,
                    "angle": 0.0,
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
                    "xle": 0.875,
                    "yle": 2.25,
                    "zle": 0.15,
                    "chord": 0.2,
                    "angle": 0.0
                }
            ]
        },
        {
            "name": "Horizontal Tail",
            "angle": 0.0,
            "mirror": true,
            "numc": 8,
            "cspace": "cosine",
            "sections": [
                {
                    "xle": 2.269,
                    "yle": 0.0,
                    "zle": 0.15,
                    "chord": 0.22,
                    "angle": -2.0,
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
                    "xle": 2.314,
                    "yle": 0.66,
                    "zle": 0.15,
                    "chord": 0.175,
                    "angle": -2.0
                }
            ]
        },
        {
            "name": "Vertical Tail",
            "angle": 0.0,
            "numc": 5,
            "cspace": "cosine",
            "sections": [
                {
                    "xle": 2.239,
                    "yle": 0.0,
                    "zle": 0.15,
                    "chord": 0.25,
                    "angle": 0.0,
                    "numb": 15,
                    "bspace": "cosine"
                },
                {
                    "xle": 2.344,
                    "yle": 0.0,
                    "zle": 0.56,
                    "chord": 0.145,
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

**Typical Python Script File:**

```python
#%% Import Dependencies
from IPython.display import display_markdown
from pyvlm_files import load_package_file
from pyvlm.outputs.msh import latticeresult_to_msh
from pyvlm.outputs.prf import latticeresult_to_prf

#%% Import Geometry
jsonfilepath = 'Test_aircraft.json'
lsys = load_package_file(jsonfilepath)

#%% Display System
display_markdown(lsys)

#%% Display Results
for case in lsys.results:
    lres = lsys.results[case]
    display_markdown(lres)

#%% Mesh File Output
lres = lsys.results['Positive 1g Cruise + 15deg Side Slip']
latticeresult_to_msh(lres, r'..\results\Test_aicraft.msh')

#%% Pessure File Output
latticeresult_to_prf(lsys, r'..\results\Test_aicraft_pressures.json')

#%% 5g Trim Case
ltrm = lsys.results['Positive 5g Dive']

#%% Plot Lift Distribution
axl = ltrm.plot_trefftz_lift_distribution()

#%% Plot Y Force Distribution
axy = ltrm.plot_trefftz_yforce_distribution()

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

![](Readme.png)