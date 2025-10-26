# LH2 Simulation Examples

This directory contains example scripts demonstrating the usage of the LH2 simulation package.

## Available Examples

### basic_usage.py

Comprehensive demonstration of core functionality:

- **Geometry Module**: Tank geometry calculations
  - Volume to height conversion for horizontal cylinders
  - Cross-sectional and surface area calculations
  - Fill percentage calculations

- **Flow Module**: Fluid flow calculations
  - Choked and non-choked compressible flow
  - Vent flow rate calculations
  - Mass flow rate predictions

- **Properties Module**: Thermophysical properties
  - Liquid hydrogen density
  - Vapor pressure calculations
  - Using polynomial approximations

- **Control Module**: Control logic
  - Pressure-driven control strategies
  - Fill regime transitions (slow/fast/topping)
  - Valve and vent control

- **Integrated Example**: Combined calculation
  - Tank venting scenario
  - Mass distribution
  - Pressure relief calculations

## Running the Examples

```bash
# From the repository root
python examples/basic_usage.py
```

## Requirements

The examples use only the core lh2sim package. CoolProp is optional - if not available, polynomial approximations will be used automatically.

## Output

The examples produce formatted console output showing:
- Input parameters
- Calculated results
- Physical units
- Control decisions

All examples use SI units:
- Pressure: Pa (displayed as bar for readability)
- Temperature: K
- Length: m
- Mass: kg
- Time: s

## Extending the Examples

To create your own example:

1. Import the required modules:
```python
from lh2sim.geometry import cyl_v_to_h
from lh2sim.flow import gas_flow
from lh2sim.properties import FluidProperties
from lh2sim.control import PressureDrivenControl
```

2. Set up your scenario parameters
3. Call the appropriate functions
4. Process and display results

See `basic_usage.py` for detailed examples.
