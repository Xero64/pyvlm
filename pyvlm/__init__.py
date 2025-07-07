from .classes import LatticeOptimum as LatticeOptimum
from .classes import LatticeResult as LatticeResult
from .classes import LatticeSystem as LatticeSystem
from .classes import LatticeTrim as LatticeTrim

USE_CUPY = False

def set_cupy(use_cupy: bool = True):
    """Set whether to use cupy for calculations."""
    global USE_CUPY
    USE_CUPY = use_cupy
    if USE_CUPY:
        print('Using cupy for calculations.')
    else:
        print('Using numpy for calculations.')
