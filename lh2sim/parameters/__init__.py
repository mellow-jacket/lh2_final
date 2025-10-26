"""
LH2 Simulation Parameters Module

This module provides parameter configuration classes for LH2 transfer simulation scenarios.
Follows the design from project_instructions.md with validation and clear structure.

Key components:
- TankParameters: Tank geometry and initial conditions
- PhysicsParameters: Thermophysical constants
- TransferParameters: Transfer control parameters
- ScenarioConfig: Complete scenario configuration

References MATLAB implementations for parameter structure but uses clean Python dataclasses.
"""

from .parameters import (
    TankParameters,
    PhysicsParameters,
    TransferParameters,
    ScenarioConfig,
    create_trailer_to_dewar_scenario,
    create_pump_driven_scenario,
    create_single_tank_venting_scenario,
)

__all__ = [
    "TankParameters",
    "PhysicsParameters",
    "TransferParameters",
    "ScenarioConfig",
    "create_trailer_to_dewar_scenario",
    "create_pump_driven_scenario",
    "create_single_tank_venting_scenario",
]
