from math import pi

from .airfoil import Airfoil, airfoil_from_dat
from .bell import Bell
from .camber import NACA4, NACA6Series
from .elliptical import Elliptical
from .mass import Mass, MassCollection, masses_from_data, masses_from_json
from .spacing import (equal_spacing, full_cosine_spacing, normalise_spacing,
                      semi_cosine_spacing)
from .stability import StabilityApproximation


def elliptical_lift_force_distribution(y: list, b: float, L: float):
    return [4/pi/b*L*(1-(2*yi/b)**2)**0.5 for yi in y]

def bell_lift_force_distribution(y: list, b: float, L: float):
    return [16/3/pi/b*L*(1-(2*yi/b)**2)**1.5 for yi in y]
