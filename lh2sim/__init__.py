"""
LH2 Simulation Package

A Python implementation of liquid hydrogen transfer simulation based on
LLNL and paper MATLAB models. Simulates mass and energy balances for
LH2 transfer between tanks with thermophysical property calculations.
"""

__version__ = "0.1.0"
__author__ = "LH2 Simulation Team"

from . import properties
from . import geometry
from . import flow
from . import control

# Optional submodules: import lazily/guarded so package import doesn't fail
_optional = []

try:
    from . import simulation

    _optional.append("simulation")
except Exception:
    # simulation is optional or still under development; don't fail import
    simulation = None

try:
    from . import parameters

    _optional.append("parameters")
except Exception:
    parameters = None

try:
    from . import visualization

    _optional.append("visualization")
except Exception:
    visualization = None

__all__ = [
    "properties",
    "geometry",
    "flow",
    "control",
] + _optional
