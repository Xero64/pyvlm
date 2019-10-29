from .camber import NACA4, NACA6Series
from .spacing import normalise_spacing, equal_spacing
from .spacing import semi_cosine_spacing, semi_elliptical_spacing
from .spacing import full_cosine_spacing, full_elliptical_spacing
from .elliptical import Elliptical
from .bell import Bell
from math import pi

def elliptical_lift_distribution(y: list, b: float, L: float):
    return [4/pi/b*L*(1-(2*yi/b)**2)**0.5 for yi in y]

def bell_lift_distribution(y: list, b: float, L: float):
    return [16/3/pi/b*L*(1-(2*yi/b)**2)**1.5 for yi in y]
