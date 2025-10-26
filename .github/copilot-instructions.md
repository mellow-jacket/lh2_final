# GitHub Copilot Instructions for LH2 Simulation Project

## Project Overview

This repository contains a Python implementation of liquid hydrogen (LH2) transfer simulation, based on LLNL and paper MATLAB models. The project provides a modular, well-tested Python codebase for simulating LH2 transfer operations between tanks with detailed thermophysical property calculations.

## Required Reading Before Starting Work

**IMPORTANT**: Before starting any task, you MUST read and understand these key documentation files:

1. **`project_instructions.md`** - Comprehensive project plan with architecture, scope, validation strategy, and implementation roadmap
2. **`FINAL_STATUS.md`** - Current project completion status, deliverables, metrics, and quality assessments
3. **`IMPLEMENTATION_SUMMARY.md`** - Technical implementation details for all modules
4. **`PROGRESS.md`** - Development progress tracker with completed/in-progress items

These files provide critical context about:
- Project goals and acceptance criteria
- Architecture decisions and module boundaries
- Physics scope and extensibility requirements
- Testing strategy and quality standards
- Current implementation status
- What's completed vs what's still needed

## MATLAB Reference Codebases

The repository includes reference MATLAB implementations for inspiration:

### MATLAB Code Locations
- **`matlab_codebases/LLNL_model/`** - Original LLNL MATLAB implementation (NASA Open Source Agreement v1.3)
- **`matlab_codebases/paper_model/`** - Enhanced derivative MATLAB code with pump-driven support

### High-Level Summaries
- **`matlab_codebases/llnl_rundown.md`** - Detailed analysis of LLNL MATLAB codebase, file-by-file breakdown
- **`matlab_codebases/paper_rundown.md`** - Analysis of paper model differences and improvements

**Usage Guidelines for MATLAB References**:
- Use for understanding physics implementation and control logic
- Use for validation targets and test case ideas
- Reference for thermophysical property correlations
- **DO NOT copy code directly** - implement clean-room Python equivalents
- Be aware of licensing (NASA OSA 1.3 for LLNL, GPL-3.0 for paper model)

## Project Architecture

### Package Structure
```
lh2sim/
├── geometry/           # Tank geometry calculations (horizontal/vertical cylinders)
├── flow/              # Fluid flow (choked/non-choked compressible flow)
├── properties/        # Thermophysical properties (CoolProp/polynomial backends)
├── control/           # Control logic (pressure-driven/pump-driven)
├── simulation/        # ODE integration (WIP)
├── parameters/        # Scenario configuration (WIP)
└── visualization/     # Plotting utilities (WIP)
```

### Key Design Principles
1. **Modular Architecture**: Clear module boundaries with explicit interfaces
2. **No Global State**: Pure functions with explicit I/O (unlike MATLAB `evalin`/`assignin`)
3. **Property Abstraction**: Backend-agnostic (CoolProp/REFPROP/polynomial)
4. **Comprehensive Testing**: Unit tests for all physics calculations
5. **Type Safety**: Type hints where beneficial
6. **Documentation**: NumPy-style docstrings throughout

### Dependencies
- **NumPy** >= 1.20.0 - Array operations
- **SciPy** >= 1.7.0 - Scientific computing, ODE integration
- **matplotlib** >= 3.3.0 - Visualization
- **CoolProp** >= 6.4.0 - Thermophysical properties (optional but recommended)
- **Pint** - Unit handling (planned)

## Code Style and Quality Standards

### Code Quality Requirements
- **Type Hints**: Use where beneficial for clarity
- **Docstrings**: NumPy-style for all public functions and classes
- **Testing**: 
  - Target: >= 60% code coverage
  - All physics calculations must have unit tests
  - Test edge cases (empty/full tanks, boundary conditions)
  - Test realistic LH2 conditions
- **Linting**: PEP 8 compliant (use flake8/black)
- **Security**: No vulnerabilities (CodeQL verified)
- **No Hard-coded Paths**: All paths must be portable

### Testing Commands
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=lh2sim --cov-report=html

# Run specific module
pytest tests/unit/test_geometry.py
```

## Physics Implementation Notes

### Thermophysical Properties
- Use `FluidProperties` class for property lookups
- Support both CoolProp and polynomial fallback backends
- Validate property calls are within valid envelope ranges
- ParaHydrogen is the reference fluid (100% para-H₂ assumption)

### Flow Calculations
- Implement choked and non-choked compressible flow
- Use directed square root helper for flow direction preservation
- Validate pressure ratios and thermodynamic consistency

### Geometry
- Support horizontal and vertical cylinders
- Volume-to-height conversions use Newton iteration
- Calculate wetted areas, interface areas, cross-sections

### Control Logic
- Implement hysteresis to prevent valve chatter
- Support both pressure-driven and pump-driven modes
- Fill regime management (slow/fast/reduced/topping)
- Vent control with deadband (±5-10% thresholds)

## Common Tasks and Guidelines

### Adding New Features
1. Check if similar functionality exists in MATLAB codebases
2. Read relevant sections of `project_instructions.md` for context
3. Design clean Python API (no global state)
4. Write unit tests first (TDD approach)
5. Implement feature with proper type hints and docstrings
6. Verify tests pass and coverage is adequate
7. Update relevant documentation files

### Modifying Existing Code
1. Run tests before making changes to establish baseline
2. Make minimal, surgical changes
3. Ensure tests still pass after changes
4. Update tests if behavior intentionally changed
5. Update docstrings if API changed

### Adding Tests
1. Follow existing test patterns in `tests/unit/`
2. Test edge cases (empty, full, boundary conditions)
3. Test realistic LH2 conditions (P=1-10 bar, T=15-30K)
4. Test physical validity (conservation, pressure ranges)
5. Use pytest fixtures for common setup

### Debugging Issues
1. Check `PROGRESS.md` for known issues
2. Reference MATLAB implementations for expected behavior
3. Verify property backend is working correctly
4. Check units and dimensional consistency
5. Validate against physics fundamentals

## Documentation Update Requirements

**CRITICAL**: Upon completion of any significant task, you MUST update these three files:

### 1. FINAL_STATUS.md
Update with:
- New deliverables completed
- Updated metrics (lines of code, test count, coverage)
- New capabilities added
- Quality metrics (security, linting status)
- Current project status summary

### 2. IMPLEMENTATION_SUMMARY.md
Update with:
- Technical details of what was accomplished
- New modules or functions added
- Implementation decisions made
- Dependencies added or changed
- Usage examples for new functionality
- Any technical debt or future work identified

### 3. PROGRESS.md
Update with:
- Mark completed items as [x]
- Add new work items as needed
- Update "In Progress" section
- Note any blockers or issues encountered
- Update testing status and coverage

**Format**: Keep consistent with existing structure using markdown checklists:
- `[x]` for completed items
- `[ ]` for pending items
- Use clear, descriptive item names

## Known Limitations and Future Work

### Completed (from PROGRESS.md)
- ✅ Geometry module (87% coverage)
- ✅ Flow module (96% coverage)
- ✅ Properties module (18% coverage - limited due to CoolProp dependency)
- ✅ Control module (95% coverage)

### High Priority (Not Yet Implemented)
- Parameters module for scenario configuration
- Simulation module for ODE integration
- Integration tests for end-to-end validation
- Event detection for vent switching

### Medium Priority
- Visualization module for plotting
- More comprehensive examples
- Performance optimization

### Low Priority
- CLI interface
- CI/CD pipeline
- Advanced plotting features

## License and IP Considerations

### Original MATLAB Code
- **LLNL model**: NASA Open Source Agreement (NOSA) v1.3
- **Paper model**: GPL-3.0

### Python Implementation
- Clean-room implementation (no code copying from MATLAB)
- Recommend BSD-3 or MIT license for Python code
- Maintain clear attribution to original work
- Keep license boundaries clear

### What You Can Do
- Read MATLAB code for understanding
- Replicate physics behavior and algorithms
- Use as validation reference
- Learn from architecture patterns

### What You Cannot Do
- Copy MATLAB code directly
- Copy comments verbatim
- Incorporate NOSA/GPL code without proper licensing
- Claim original authorship of physics models

## Special Notes

### Property Backend Considerations
- Code should work with or without CoolProp installed
- Use polynomial fallbacks for testing/CI
- CoolProp recommended for production accuracy
- Support for REFPROP may be added in future (optional, licensed)

### MATLAB Translation Pitfalls
- **Global State**: MATLAB code uses `evalin('base',...)` - eliminate this pattern
- **Duplicate Functions**: MATLAB has duplicated helpers - consolidate in Python
- **Units**: MATLAB mixes units - use consistent SI in Python (consider Pint)
- **REFPROP Dependency**: Abstract behind property backend interface

### Performance Considerations
- Property evaluations are hotspots
- Event detection can be expensive
- Consider vectorization where possible
- Use appropriate ODE solver (scipy.integrate.solve_ivp)
- Process-based parallelism for parameter sweeps

## Quick Reference Commands

```bash
# Setup
pip install -r requirements.txt
pip install -r requirements-dev.txt
pip install -e .

# Testing
pytest
pytest --cov=lh2sim --cov-report=html
pytest tests/unit/test_geometry.py -v

# Linting
flake8 lh2sim/
black lh2sim/

# Run examples
python examples/basic_usage.py
```

## Getting Help

1. Read `project_instructions.md` for comprehensive project context
2. Check `PROGRESS.md` for current status and known issues
3. Review MATLAB rundown files for physics understanding
4. Examine existing tests for usage patterns
5. Look at completed modules for code style examples

## Summary

This project is building a modern Python LH2 simulation framework that:
- Replicates MATLAB model capabilities with clean Python design
- Provides modular, testable, well-documented code
- Supports multiple thermophysical property backends
- Enables reproducible scientific simulations
- Maintains clear IP boundaries and licensing

Always read the key documentation files first, update them upon completion, and refer to MATLAB codebases for inspiration (not copying). Maintain high code quality standards and comprehensive testing throughout.
