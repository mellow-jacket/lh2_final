"""
LH2 Simulation Module

This module provides the core simulation engine for LH2 transfer operations.
It implements mass and energy balance equations, event detection, and ODE integration.

Key components:
- SimulationState: State vector management
- Simulator: Main simulation orchestrator
- mass_balance: Mass balance equations
- energy_balance: Energy balance equations (stub for future)

This is a simplified initial implementation focusing on mass balance.
Energy balance and detailed heat transfer will be added in future iterations.
"""

from .simulation import (
    SimulationState,
    SimulationResult,
    Simulator,
)

__all__ = [
    "SimulationState",
    "SimulationResult",
    "Simulator",
]
