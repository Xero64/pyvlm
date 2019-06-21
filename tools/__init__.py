from math import pi

def elliptical_lift_distribution(y: list, b: float, L: float):
    return [4/pi/b*L*(1-(2*yi/b)**2)**0.5 for yi in y]

def bell_lift_distribution(y: list, b: float, L: float):
    return [16/3/pi/b*L*(1-(2*yi/b)**2)**1.5 for yi in y]
