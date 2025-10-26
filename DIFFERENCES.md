# DIFFERENCES.md — refreshed after recent updates

This document records the current differences between the LLNL / paper MATLAB reference models and the Python port in `lh2sim/`, reflecting the recent updates you made (added energy-balance helpers, expanded simulation scaffolding, partial property updates, etc.). It lists what is implemented, what's partially implemented, and what's missing — with prioritized next steps.

Summary of today's observation
- The repository has advanced: `lh2sim` now contains more complete scaffolds and several concrete modules (including a new `energy_balance` module under `lh2sim/simulation`).
- The Python code is moving from a skeleton to a functional scientific codebase, but the remaining work is primarily scientific correctness (thermophysical property parity, internal-energy-based energy balances, and control/event integration).

Detailed status by subsystem

1) Properties (thermophysical)
- Current state: partial / improved
  - `lh2sim/properties/properties.py` implements a `FluidProperties` abstraction, CoolProp detection, and a `vapor_pressure` function that references the original MATLAB polynomial approach. The file includes placeholders and comments indicating polynomial coefficients are used.
  - Action required: verify polynomial coefficients, confirm unit conventions (Pa vs kPa), and add automated tests comparing the polynomial fallback and CoolProp against a vetted reference (REFPROP if available, otherwise a saved reference table).

2) Simulation & energy balance
- Current state: scaffolded + helper module added
  - `lh2sim/simulation/simulation.py` now contains `Simulator` with `_compute_initial_state`, `_derivatives`, `_enrich_result`, `_create_event_functions`, and `run_with_events` methods. The code imports `energy_balance` helpers.
  - `lh2sim/simulation/energy_balance.py` exists and provides routines to assemble multi-node grids, compute latent heat, and helper heat-transfer correlations (these functions are partially implemented and need to be integrated into `_derivatives`).
  - Action required: wire energy_balance functions into `_derivatives`, implement node-by-node internal-energy equations (MATLAB uses U as state), and ensure mass/energy conservation.

3) Geometry
- Current state: implemented
  - `lh2sim/geometry/geometry.py` implements `cyl_v_to_h` and geometric helpers. Add unit tests for edge cases and numeric parity with MATLAB.

4) Flow
- Current state: implemented/partial
  - `lh2sim/flow/flow.py` implements `gas_flow`, `directed_sqrt`, `valve_flow_coefficient`, and `vent_flow_rate`. Structure mirrors MATLAB `gasFlow.m`.
  - Action required: numeric unit tests for choked vs non-choked regime, and sign-handling tests for flow reversal.

5) Control
- Current state: partial/implemented
  - `lh2sim/control/control.py` includes `PressureDrivenControl` and `PumpDrivenControl` classes with hysteresis state and interface `compute_control`. Some helper methods exist but bodies are omitted in places.
  - Action required: port exact setpoints/deadbands from MATLAB `Parameters_*.m` and finish controller logic (esp. vent & vaporizer decision trees and pump-mode behavior).

6) Parameters & scenarios
- Current state: partial
  - `lh2sim/parameters/parameters.py` defines dataclasses and scenario factory functions (trailer→dewar, pump-driven) with several default values. Not all MATLAB `Parameters_*.m` files are ported.
  - Action required: port remaining parameter files and unify naming/units.

7) Visualization & output
- Current state: partial
  - `lh2sim/visualization/visualization.py` provides plotting scaffolds for tank levels, pressures, temperatures, masses, densities, and a dashboard. Detailed MATLAB plot annotations and the exact node-level export format remain to be matched.
  - Action required: implement export routines and ensure plotted series map to the same state variables/indices as the MATLAB output.

8) I/O, logging, and utilities
- Current state: missing/low priority
  - No direct equivalents for MATLAB's waitbar, `UpdateXLSLog`, or text export orchestration; recommend adding `lh2sim/io` or `lh2sim/utils` for CSV/Excel export and long-run progress reporting.

9) Tests and validation
- Current state: incomplete/required
  - Add unit tests for `vapor_pressure`, `cyl_v_to_h`, and `gas_flow`. Add an integration smoke test running `Simulator` for a short time horizon and asserting mass and energy conservation within tolerances.

10) Licensing & provenance
- Current state: informational
  - MATLAB LLNL model is NASA OSA 1.3. Confirm README states the MATLAB code is a reference and that Python is clean-room unless any text/code was copied with compatible terms.

Prioritized next steps (short list — pick one to start)
1. Verify and unit-test `vapor_pressure` polynomial coefficients and add reference table for CI (recommended first).
2. Implement node-by-node energy balances in `_derivatives` using `energy_balance` helpers and run a smoke integration test.
3. Port remaining `Parameters_*.m` values into `ScenarioConfig` and wire to controllers; add hysteresis unit tests.
4. Add unit tests for `gas_flow` (choked/non-choked and reversal) and `cyl_v_to_h` edge cases.
5. Implement results export routine that reproduces MATLAB node-level CSV/TXT layout and finish plotting parity.

please complete items in the next steps
