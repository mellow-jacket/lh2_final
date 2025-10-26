"""
Properties Module

Provides thermophysical property calculations for liquid hydrogen using CoolProp.
Includes abstraction layer to support both CoolProp and REFPROP backends.

Key functionality:
- Density, enthalpy, entropy calculations
- Vapor pressure and saturation properties
- Transport properties (viscosity, thermal conductivity)
- Polynomial correlations as fallbacks
"""

from .properties import (
    FluidProperties,
    vapor_pressure,
    COOLPROP_AVAILABLE,
)

__all__ = [
    "FluidProperties",
    "vapor_pressure",
    "COOLPROP_AVAILABLE",
]
