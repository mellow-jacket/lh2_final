# Differences and Missing Features

This document catalogs differences and missing features between the original LLNL MATLAB codebase (and the derived "paper" codebase) and the current Python implementation in `lh2sim/` found in this repository. It is intended as a checklist and guide for prioritizing implementation, testing, and validation work needed to reach feature parity.

NOTE: this report was generated from the MATLAB files provided (LLNL_model and paper_model summaries) and the Python modules under `lh2sim/` (simulation, properties, geometry, flow, control, parameters, visualization). Where the Python code contains large omitted blocks ("/* Lines omitted */") those are treated as not implemented unless clearly present.

## High-level summary

- The Python port contains a good modular skeleton: `simulation`, `properties`, `geometry`, `flow`, `control`, `parameters`, and `visualization` modules.
- The MATLAB LLNL and paper models implement a mature, multi-node (layered) thermal model (liquid and vapor nodes per tank), REFPROP usage via `refpropm`, detailed energy balances, many parameter sets, plotting/saving utilities, and scripts that orchestrate full scenarios.
- The Python port currently implements mass-balance-focused ODE scaffolding and property-call abstractions, but several physics, numerical, I/O, and control behaviors from the MATLAB code are missing or incomplete. The sections below list these differences in detail and map them to files where possible.

## Methodology & assumptions

- I compared the MATLAB LLNL files you provided (notably `LH2Simulate.m`, `LH2Control.m`, `vaporpressure.m`, `cylVToH.m`, `plotLH2Data.m`, `runNominal.m`, `inputs_*`, and `Parameters_*.m`) with the Python sources in `lh2sim/` (where many functions are partially stubbed or have omitted blocks).
- Where the Python file contains large omitted blocks, I treated those as missing functionality.
- I assume the MATLAB code uses REFPROP (`refpropm`) and extensive array-based state vectors (mL1, mV1, mL2, mV2, layer temperatures, wall nodes, etc.). The Python port uses a simplified state and will require expansion to match arrays and energy balances.

## Detailed differences and missing features

Below each section lists the missing or different behavior, followed by file pointers (MATLAB → Python) and suggested priority.

### 1) Thermophysical property backend parity

- MATLAB: LLNL and paper models call `refpropm` (REFPROP) extensively for U, P, D, h, etc., with units handling and parahydrogen options.
- Python: `lh2sim/properties.py` has a `FluidProperties` abstraction with CoolProp backend and polynomial fallbacks. However:
  - REFPROP support is not implemented/verified. The MATLAB uses REFPROP-specific behavior and units; parity requires either REFPROP bindings (if licensed) or validated CoolProp/empirical polynomials with test coverage.
  - The polynomial fallback functions appear minimal and placeholders in Python (e.g., `_liquid_density_polynomial`, `_vapor_pressure_polynomial` are stubbed or simplified). The MATLAB `vaporpressure.m` contains specific polynomial coefficients and a branch for supercritical/scenario — that specific logic is only partially reproduced in `properties.vapor_pressure`.
  - Many property call signatures used in MATLAB (U from T and Q, P from T and D, etc.) need exact equivalents with matching units and numeric tolerances.

Priority: high — property accuracy underpins the entire physics model and numerical stability.

Files: MATLAB `vaporpressure.m`, `LH2Simulate.m` → Python `lh2sim/properties/properties.py`.

### 2) Multi-node (layered) mass & energy balances and stratification

- MATLAB: Tanks are discretized into arrays of liquid nodes (`nL1`, `nL2`) and vapor nodes (`nV1`, `nV2`), each node with its own temperature and internal energy (`UL10`, `Uv10`, etc.). Energy balances and heat flux between nodes (and with walls) are present.
- Python: `simulation.Simulator` and `SimulationState` are currently simplified; their dataclasses are empty or incomplete and the `_derivatives` implementation is mostly omitted. The following are missing or incomplete:
  - Node arrays for multi-cell liquid and vapor temperatures and energies.
  - Heat transfer between layers, wall nodes, and boil-off/condensation models.
  - Liquid/vapor phase tracking per-node and consistent phase-change mass transfer terms.
  - Proper initial-condition construction that mirrors MATLAB's `x0` vector with mass and internal energy per node.

Priority: high — needed to replicate stratification, boil-off, and accurate thermal dynamics.

Files: MATLAB `LH2Simulate.m` (initialization and state vector) → Python `lh2sim/simulation/simulation.py`.

### 3) Energy balances and full thermodynamic state

- MATLAB: Explicit energy balances via REFPROP internal energy (`U`) and enthalpy calls used in mass/energy flux computations. Terms like Qdot, JL2*(refprop('U', ...)+P.refstatedelta) appear.
- Python: Energy-balance functions are stubbed — simulation currently focuses on mass balances. The `energy_balance` is a placeholder. Implementation details such as work terms, heat leaks, condensation/evaporation latent heat are missing.

Priority: high.

Files: MATLAB `LH2Simulate.m`, `vaporpressure.m`, energy-related lines → Python `lh2sim/simulation` and `lh2sim/properties`.

### 4) ODE state vector layout, event detection & integration parity

- MATLAB: creates a single state vector `x0` with many components (masses, energies, possibly pressures and temperatures), uses MATLAB ODE solver with event detection for fill completion and vent switching, and uses global variables for control flags.
- Python: `Simulator.run` uses `scipy.integrate.solve_ivp` but `_compute_initial_state` and `_derivatives` are incomplete. Event functions and stopping conditions (e.g., ET_fill_complete, ST_vent_complete) are not implemented.

Missing items:
  - Exact state vector mapping to MATLAB `x0` (indices & units).
  - Event functions that set flags and stop integration when required.
  - Reproducible control interface to call controller logic each timestep with the same hysteresis behavior.

Priority: high.

Files: MATLAB `LH2Simulate.m` (event detection) → Python `lh2sim/simulation/simulation.py`.

### 5) Control logic differences (venting, vaporizer, fill regimes, hysteresis)

- MATLAB: `LH2Control.m` implements regime-dependent logic based on ET liquid height, ST/ET pressures, and uses `getSTVentState` with hysteresis (threshold ±5% deadbands). The paper model includes both pressure-driven and pump-driven versions with different thresholds and pump mass flow rates.
- Python: `lh2sim/control/control.py` contains `PressureDrivenControl` and `PumpDrivenControl`, but many internals are omitted (hysteresis threshold values, mapping of ET vs ST fill regimes, the exact vent setpoints `p_ST_slow`, `p_ST_fast`, `p_ST_final`, and the ET fill completion logic). Pump behavior (pump speed control, startup blowdown, plug power) from paper model is not fully modeled.

Missing items:
  - Exact setpoints and deadbands from `Parameters_*.m`.
  - Implementation of `ET_fill_complete` and `ST_vent_complete` flags consistent with MATLAB.
  - Vaporizer valve logic implemented exactly (open fraction, hysteresis, PID-like behavior if present).

Priority: high/medium depending on use-case (pressure-driven scenarios higher).

Files: MATLAB `LH2Control.m`, `LH2Control_pump.m`, parameters files → Python `lh2sim/control/control.py` and `lh2sim/parameters/parameters.py`.

### 6) Flow and orifice models: choked flow and careful sign handling

- MATLAB: `gasFlow.m` handles choked/non-choked compressible flow with precise isentropic relations and signed square-root helpers.
- Python: `flow.gas_flow` implements choked/non-choked logic but many details are omitted. Verify constants, critical pressure ratio expression, and sign/direction conventions. `directed_sqrt` and valve area routines exist but must be validated against MATLAB formulas and units.

Missing items:
  - Exact constants and downstream vs upstream density/temperature choices used in MATLAB.
  - Consistent handling when P1 < P2 (flow reversal) for mass flow sign and upstream property selection.

Priority: medium.

Files: MATLAB `gasFlow.m` → Python `lh2sim/flow/flow.py`.

### 7) Geometry: cylVToH differences and numerical method details

- MATLAB: `cylVToH.m` is an implementation with Newton iteration and particular initial guesses and branches for s > A/2. The code handles near-zero and near-full cases, uses `atan` expressions and specific tolerances.
- Python: `geometry.cyl_v_to_h` replicates the approach, but the omitted lines may hide differences in initial guess, variable mapping (x vs height), and edge-case handling (exact half-full check, small volume behavior). Ensure convergence tolerance and mapping to H (R ± x) consistent with MATLAB.

Priority: low/medium — geometry likely fine but verify numeric parity.

Files: MATLAB `cylVToH.m` → Python `lh2sim/geometry/geometry.py`.

### 8) Vapor pressure and phase determination specifics

- MATLAB `vaporpressure.m` contains a polynomial mapping from vapor density to saturation temperature, branch for supercritical state (calls REFPROP for single-phase) and a polynomial expression for two-phase pv. The coefficients are explicit.
- Python `properties.vapor_pressure` contains a simplified placeholder (T_sat = 20.0 + 0.1*rho_vapor) and comments that the original coefficients are needed.

Missing items:
  - Exact polynomial coefficients and logic from MATLAB to replicate phase boundary and consistent pv calculation.
  - Unit conventions: MATLAB returns Pa (with multiplication by 1e3 in places). The Python function must match units.

Priority: high.

Files: MATLAB `vaporpressure.m` → Python `lh2sim/properties/properties.py`.

### 9) Output, plotting, and result extraction parity

- MATLAB: `plotLH2Data.m` produces many subplots, predefined markers (70%, 80%, 90% lines), and saves data with `dlmwrite`. Scripts `runNominal.m` and `MAIN.m` orchestrate running, extracting data, plotting, and saving.
- Python: `visualization.visualization.py` provides plotting functions but many lines are omitted and likely do not implement the detailed dashboard and all figures from MATLAB. Also saving results in MATLAB includes writing many arrays (TL2, Ts2, Tv2, Tw2, uL2, uv2, etc.) that come from multi-node state arrays.

Missing items:
  - Parity of plots and exported CSV/txt including all node-level temperature/time series, boiloff metrics, and energy budgets.
  - Recreated 'extract data' step that produces the numerous columns MATLAB writes to disk.

Priority: medium.

Files: MATLAB `plotLH2Data.m`, `runNominal.m` → Python `lh2sim/visualization/visualization.py` and top-level examples in `examples/`.

### 10) Parameters and scenario sets

- MATLAB: Multiple `Parameters_*.m` files define a wide range of constants and scenario-specific presets (press-diff vs pump, initial blowdown, plug power, geometric details in gallons or liters to m³ conversion, parahydrogen settings, p_ST setpoints, pump flow rates, Topfill variables, etc.).
- Python: `parameters.parameters.py` includes dataclasses and two scenario factory functions, but many specific parameter values from MATLAB are omitted (marked with "Lines omitted").

Missing items:
  - Complete translation of `Parameters_*.m` values into Python `ScenarioConfig` and per-method presets.
  - Mapping of MATLAB naming conventions (P.* fields) to the Python dataclasses.

Priority: medium/high — to reproduce behavior for validation scenarios.

Files: MATLAB `Parameters_*.m` → Python `lh2sim/parameters/parameters.py`.

### 11) Pump-driven variant behavior (paper model)

- MATLAB paper model includes pump-driven flows (`LH2Simulate_Pump.m`, `LH2Control_pump.m`, special parameter files). The paper model implements pump mass transfer rates, possible initial blowdown logic, plug power, and a slightly different control strategy.
- Python: `control.PumpDrivenControl` exists but is not fully implemented (internal logic omitted), and parameter mapping for pump rates and pump-specific event handling is incomplete.

Priority: medium (if pump-mode scenarios are required for your work, raise to high).

Files: MATLAB `LH2Simulate_Pump.m`, `LH2Control_pump.m`, `Parameters_TrailerToMain_Pump*.m` → Python `lh2sim/control` + `lh2sim/parameters`.

### 12) I/O, logging, waitbars, and interactive scripts

- MATLAB: `waitbar`, `display`, `dlmwrite`, `save` and other helpful scripts exist for user interaction and data export. The paper model includes `waitbartime`, logging to Excel via `UpdateXLSLog.m`, and save routines.
- Python: There is no direct equivalent for interactive waitbars or `UpdateXLSLog`. `examples/` exist but exporting and logging utilities are not at MATLAB parity.

Priority: low.

Files: MATLAB `waitbartime.m`, `UpdateXLSLog.m` → consider `lh2sim/utils` or `examples/` enhancements.

### 13) Testing and validation harness

- MATLAB: implicit validation via reference runs and derivation from REFPROP/LLNL code. The repository includes many tests for the Python side (per workspace tree, `tests/unit/` shows some tests), but coverage may be incomplete for properties and energy/flow edge cases.
- Python: `tests/` exists but property-level unit tests for polynomial coefficients, `vapor_pressure`, and multi-node energy balance are missing or limited.

Priority: high — property and energy tests are required for scientific correctness.

Files: `tests/unit/*` ↔ new tests to add for `properties`, `geometry`, `flow`, `simulation` event detection, and `control` hysteresis.

### 14) Licensing & metadata

- MATLAB LLNL model is NASA OSA 1.3 (contained in `LLNL_model/LICENSE.md`). The Python repo should maintain a clear license and note that MATLAB files are reference only and not copied. Ensure license compatibility if incorporating any parts of MATLAB code.

Priority: informational/legal — ensure compliance.

Files: `matlab_codebases/LLNL_model/LICENSE.md` → update repo docs if needed.

## Concrete prioritized next steps (recommended)

1. Property backend and vapor pressure parity (High)
   - Add/refine polynomial coefficients from `vaporpressure.m` and validate against REFPROP or CoolProp for a range of states. Implement REFPROP backend optionally (if license available) or document divergence.
   - Add unit tests for property functions covering the range used in scenarios.

2. Expand state vector & energy balances (High)
   - Implement multi-node liquid/vapor arrays, wall nodes, heat-leak terms, and energy derivative calculations mirroring MATLAB `x0` and `_derivatives` logic.
   - Create tests that replicate a known MATLAB run (or nominal case) to compare time series.

3. Implement control logic parity (High/Medium)
   - Translate `LH2Control.m` and `LH2Control_pump.m` hysteresis and thresholds exactly. Add parameter sets from `Parameters_*.m`.

4. Implement event detection and flags (High)
   - Reproduce `ET_fill_complete`, `ST_vent_complete`, and other event-based behaviors. Hook them to `solve_ivp` events.

5. Validate `gas_flow` and flow reversal handling (Medium)
   - Add unit tests for choked vs non-choked and ensure sign/direction matches MATLAB.

6. Finish visualization parity and data export (Medium)
   - Add full dashboard plots and an export routine that reproduces the MATLAB `dlmwrite` output layout (node-level timeseries).

7. Tests and CI (High)
   - Create focused unit tests and an integration test that runs a short scenario and compares a few key metrics to MATLAB reference (or known expected values).

8. Documentation and examples (Low-Medium)
   - Update README and add a migration doc describing differences in physics assumptions (CoolProp vs REFPROP).

## Appendix: File mapping (quick reference)

- MATLAB LLNL/paper → Python
  - `LH2Simulate.m`, `LH2Simulate_Pump.m` → `lh2sim/simulation/simulation.py`
  - `LH2Control.m`, `LH2Control_pump.m` → `lh2sim/control/control.py`
  - `vaporpressure.m` → `lh2sim/properties/properties.py` (function `vapor_pressure`)
  - `cylVToH.m` → `lh2sim/geometry/geometry.py` (function `cyl_v_to_h`)
  - `gasFlow.m` → `lh2sim/flow/flow.py` (function `gas_flow`)
  - `plotLH2Data.m` → `lh2sim/visualization/visualization.py`
  - `Parameters_*.m` → `lh2sim/parameters/parameters.py` (dataclasses & scenario factories)

## Closing notes

I created this high-detail checklist to guide targeted implementation work. If you want, I can:

- (A) Start implementing a high-priority item (property polynomial coefficients and unit tests) and validate against CoolProp.
- (B) Implement the multi-node initial state and a minimal energy-balance term to get an end-to-end example running and produce plots comparable to MATLAB for a short time horizon.

Tell me which of the prioritized items you'd like me to implement first and I'll begin with a focused plan and corresponding PR-style patch.

---
Generated on: 2025-10-26
