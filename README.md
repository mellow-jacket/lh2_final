# LH2 Simulation - Python Implementation

A Python implementation of liquid hydrogen (LH2) transfer simulation based on LLNL and paper MATLAB models. Simulates mass and energy balances for LH2 transfer between tanks with detailed thermophysical property calculations.

## Overview

This project provides a modular, well-tested Python codebase for simulating liquid hydrogen transfer operations. It includes:

- **Thermophysical Properties**: CoolProp integration with polynomial fallbacks
- **Flow Calculations**: Choked/non-choked compressible flow models
- **Geometry**: Tank geometry calculations for horizontal and vertical cylinders
- **Control Logic**: Pressure-driven and pump-driven control strategies
- **Simulation**: ODE-based time integration with event handling
- **Visualization**: Comprehensive plotting utilities

## Installation

### Requirements
- Python >= 3.8
- NumPy >= 1.20.0
- SciPy >= 1.7.0
- CoolProp >= 6.4.0
- matplotlib >= 3.3.0

### Install from source

```bash
# Install dependencies
pip install -r requirements.txt

# Install development dependencies (for testing)
pip install -r requirements-dev.txt

# Install package in development mode
pip install -e .
```

## Quick Start

```python
import numpy as np
from lh2sim.geometry import cyl_v_to_h
from lh2sim.flow import gas_flow
from lh2sim.properties import FluidProperties

# Calculate liquid height in horizontal cylinder
R = 1.5  # m, cylinder radius
L = 10.0  # m, cylinder length
V = 20.0  # m³, liquid volume
height = cyl_v_to_h(V, R, L)
print(f"Liquid height: {height:.2f} m")

# Calculate gas flow through orifice
CA = 0.001  # m², effective flow area
rho = 1.3  # kg/m³, gas density
P1 = 1.5e5  # Pa, upstream pressure
P2 = 1.0e5  # Pa, downstream pressure
gamma = 1.4  # ratio of specific heats
mdot = gas_flow(CA, rho, P1, P2, gamma)
print(f"Mass flow rate: {mdot:.4f} kg/s")

# Get thermophysical properties
props = FluidProperties("Hydrogen", backend="CoolProp")
T = 20.0  # K
P = 1e5  # Pa
rho_liquid = props.density(T=T, P=P)
print(f"Liquid density at {T}K: {rho_liquid:.2f} kg/m³")
```

## Testing

Run the test suite:

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=lh2sim --cov-report=html

# Run specific test module
pytest tests/unit/test_geometry.py
```

## Project Structure

```
lh2_final/
├── lh2sim/                 # Main package
│   ├── geometry/           # Geometric calculations
│   ├── flow/              # Flow calculations
│   ├── properties/        # Thermophysical properties
│   ├── control/           # Control logic (WIP)
│   ├── simulation/        # ODE simulation (WIP)
│   ├── parameters/        # Scenario parameters (WIP)
│   └── visualization/     # Plotting utilities (WIP)
├── tests/                 # Test suite
│   └── unit/             # Unit tests
├── LLNL_model/           # Original MATLAB code (reference)
├── paper_model/          # Derived MATLAB code (reference)
├── PROGRESS.md           # Development progress tracker
└── pyproject.toml        # Package configuration
```

## Status

**Current Version**: 0.1.0 (In Development)

See `PROGRESS.md` for detailed development status and roadmap.

### Completed
- ✅ Geometry module with unit tests
- ✅ Flow module with unit tests
- ✅ Properties module with CoolProp integration
- ✅ Comprehensive test suite

### In Progress
- 🔄 Control module
- 🔄 Simulation module
- 🔄 Parameters module
- 🔄 Visualization module

## Background

Based on MATLAB codebases from LLNL and paper models. See `llnl_rundown.md` and `paper_rundown.md` for detailed analysis.

## License

See LICENSE file. Original MATLAB code under NASA Open Source Agreement (NOSA) v1.3.


# generate phase diagram
python -m lh2sim.properties.plot_phase_diagram