# Differences and Missing Features (refreshed)

This document catalogs differences and missing features between the original LLNL MATLAB codebase (and the derived "paper" codebase) and the current Python implementation in `lh2sim/` found in this repository. It is a refreshed, up-to-date comparison produced after the latest Python updates. Use it to prioritize follow-up work.

Note about source material: I compared the provided MATLAB LLNL and paper_model attachments (key files: `LH2Simulate.m`, `LH2Control.m`, `vaporpressure.m`, `cylVToH.m`, `gasFlow.m`, `plotLH2Data.m`, `Parameters_*.m`, etc.) with the current Python modules under `lh2sim/`. The Python files contain a mix of concrete implementations and omitted/stubbed sections; where logic is present I mark it "implemented/partial" and where it's missing I mark it "missing/needs work."

## High-level summary (updated)

- The Python package `lh2sim` now exposes `properties`, `geometry`, `flow`, and `control` and conditionally imports `simulation`, `parameters`, and `visualization` (guarded in `lh2sim/__init__.py`).
- Several core algorithms have scaffolding or partial implementations: `geometry.cyl_v_to_h`, `flow.gas_flow` and helpers, `properties.vapor_pressure` (now references MATLAB polynomial approach), `control` classes for pressure/pump-driven strategies, and `simulation.Simulator` contains event-related method shells.
- Remaining high-impact gaps are scientific: full multi-node energy equations, complete property parity with REFPROP, event and controller integration parity, and a rigorous test harness to validate numerics against reference data.

## Detailed differences and current status

I list each major area, what's now implemented/partial, and what's still missing.

1) Thermophysical property backend parity
- Status: partial/implemented
  - `lh2sim/properties/properties.py` has a `FluidProperties` abstraction and a `vapor_pressure` function referencing the MATLAB polynomial approach. The module detects CoolProp and will use it if available.
  - Missing/To-verify: explicit REFPROP parity is still absent (REFPROP requires license). The polynomial fallbacks need explicit coefficient verification, and unit tests should be added to validate the polynomials and CoolProp results across the operating envelope.

2) Multi-node mass & energy balances and stratification
- Status: scaffolded but incomplete
  - `lh2sim/simulation/simulation.py` contains `Simulator`, `_compute_initial_state`, `_derivatives`, `_enrich_result`, `_create_event_functions`, and `run_with_events` placeholders — the solver integration (solve_ivp) is present.
  - Missing: the concrete multi-cell liquid/vapor node arrays, energy-balance algebra using internal energy (U), latent heat treatment, inter-node heat transfer, and proper initial-state mapping of MATLAB's `x0` vector.

3) Energy balances and full thermodynamic state
- Status: missing/partial
  - Energy balance equations that rely on node-specific internal energies are not yet implemented. The `simulation` module currently focuses on mass-balance scaffolding.
  - Implementing internal-energy-based ODE terms and latent heat terms is required for correct thermal dynamics.

4) ODE state-vector layout, event detection & integration parity
- Status: partial
  - Event scaffolding exists (`_create_event_functions`, `run_with_events`), and `solve_ivp` usage is in place.
  - Missing: full event functions (fill complete, vent open/close transitions), mapping the event outcomes to controller flags, and a documented state-index mapping consistent with MATLAB.

5) Control logic (venting, vaporizer, fill regimes, hysteresis)
- Status: partial/implemented
  - `PressureDrivenControl` and `PumpDrivenControl` classes exist with internal hysteresis state and helper methods. These offer a clean interface (`compute_control`) to call from the simulator.
  - Missing: exact setpoints and deadbands from MATLAB `Parameters_*.m` and some internal control decision logic that was in MATLAB (some methods in the Python controller have omitted bodies). Pump-mode specifics (blowdown, plug power) require porting.

6) Flow and orifice models: choked flow and sign handling
- Status: implemented/partial
  - `lh2sim/flow/flow.py` contains `gas_flow`, `directed_sqrt`, `valve_flow_coefficient`, and `vent_flow_rate` (helper). The general structure mirrors the MATLAB implementation.
  - Missing: numeric verification (unit tests) to ensure choked-flow cutoffs, coefficient usage, and sign conventions exactly match MATLAB behaviour.

7) Geometry: cylVToH and related helpers
- Status: implemented
  - `lh2sim/geometry/geometry.py` provides `cyl_v_to_h`, cross-section and lateral area functions. Newton iteration logic is implemented similar to MATLAB's `cylVToH.m`.
  - Action: add unit tests for edge cases (near-empty, half-full, near-full) to confirm numeric parity.

8) Vapor pressure and phase determination specifics
- Status: partial/implemented
  - `vapor_pressure` function in `properties.py` uses the MATLAB polynomial approach scaffold and includes handling for small/negative densities. Verify the polynomial coefficients were ported correctly and units match (Pa).

9) Output, plotting, and result extraction parity
- Status: partial
  - `lh2sim/visualization/visualization.py` includes plotting functions (tank levels, pressures, temperatures, masses, densities, summary dashboard). Many plotting scaffolds exist but detailed layout/annotation parity with MATLAB `plotLH2Data.m` is not fully verified.
  - Missing: export routine that reproduces MATLAB's node-level timeseries column layout and `runNominal`-style orchestration for saving results.

10) Parameters and scenario sets
- Status: partial
  - `lh2sim/parameters/parameters.py` contains dataclass scaffolds and factory functions for trailer-to-dewar and pump-driven scenarios with some default values present.
  - Missing: full port of all MATLAB `Parameters_*.m` files and consistent naming mapping. Some scenario parameters are still omitted.

11) Pump-driven variant behavior (paper model)
- Status: partial
  - `PumpDrivenControl` class exists; parameters for pump mass transfer are acknowledged.
  - Missing: pump ramping, initial blowdown behavior, plug-power effects, and full interaction with multi-node energy balances.

12) I/O, logging, waitbars, and interactive scripts
- Status: missing/low priority
  - No direct equivalent of MATLAB's `waitbar`, `UpdateXLSLog`, or `dlmwrite` orchestration. Consider adding `lh2sim/io` or `lh2sim/utils` helpers (CSV/Excel writing, progress reporting via tqdm).

13) Testing and validation harness
- Status: incomplete/required
  - `tests/` exists, but more numeric tests are needed: property function validation (vapor pressure, density), geometry edge-case tests, choked-flow tests, controller hysteresis tests, and at least one short integration test that asserts mass/energy conservation.

14) Licensing & metadata
- Status: informational
  - MATLAB LLNL model is NASA OSA 1.3 (see `matlab_codebases/LLNL_model/LICENSE.md`). Document the provenance and confirm the Python implementation is a clean-room implementation or note any direct derivations.

## Prioritized next steps (recommended)

1) Property backend & vapor-pressure validation (High)
   - Verify polynomial coefficients in `vapor_pressure` and add wide-range unit tests comparing to CoolProp and/or REFPROP reference values (REFPROP optional). Create a small reference dataset so CI tests are reproducible without REFPROP.

2) Implement multi-node energy balances & initial-state mapping (High)
   - Complete `SimulationState` dataclass, implement `_compute_initial_state` mapping from scenario parameters to multi-node masses/energies, and implement `_derivatives` energy terms including latent heat and inter-node conduction.

3) Complete controllers & port parameters (High)
   - Port the exact setpoints/deadbands from MATLAB `Parameters_*.m` into Python `ScenarioConfig`. Complete controller logic for venting and vaporizer behavior and add unit tests for hysteresis cases.

4) Event handling & run_with_events (High)
   - Implement event detection functions (ET_fill_complete, ST_vent_complete, etc.) and ensure `solve_ivp` events change controller state and can stop the integration as in MATLAB.

5) Flow model verification (Medium)
   - Unit tests for `gas_flow` covering choked/non-choked and reversal cases.

6) Visualization & export parity (Medium)
   - Finish plotting details and implement a `results_export` function to write node-level timeseries in the same column order as MATLAB reference runs.

7) Tests & CI (High)
   - Add focused unit tests for properties, geometry, flow, and a smoke integration test; configure CI to run these (CoolProp optional).

8) Documentation & licensing notes (Low-Medium)
   - Update `README.md` and `IMPLEMENTATION_SUMMARY.md` with current status and licensing provenance.

## Appendix: file mapping (quick reference)

- MATLAB LLNL/paper → Python (brief status)
  - `LH2Simulate.m`, `LH2Simulate_Pump.m` → `lh2sim/simulation/simulation.py` (scaffolding + event hooks present; energy equations incomplete)
  - `LH2Control.m`, `LH2Control_pump.m` → `lh2sim/control/control.py` (controller classes present; some logic pending)
  - `vaporpressure.m` → `lh2sim/properties/properties.py` (`vapor_pressure` partially ported; verify coefficients)
  - `cylVToH.m` → `lh2sim/geometry/geometry.py` (`cyl_v_to_h` implemented)
  - `gasFlow.m` → `lh2sim/flow/flow.py` (`gas_flow` implemented; verify numerics and sign)
  - `plotLH2Data.m` → `lh2sim/visualization/visualization.py` (plot scaffolds exist)
  - `Parameters_*.m` → `lh2sim/parameters/parameters.py` (dataclasses + some defaults; many parameter files remain to port)

## Closing notes

Good progress: your recent updates added important scaffolding and several concrete function implementations. The remaining high-risk work is scientific (properties + energy-balance). I recommend starting with property validation or the multi-node energy-balance implementation — both are high priority and will make subsequent controller/visualization work reliable.

Pick one of the following and I'll start immediately with a focused plan and implementation:
- Verify & test `vapor_pressure` polynomial and add unit tests (recommended first).
- Implement the multi-node energy-balance and finish `_derivatives` with a smoke-run test.
- Port all MATLAB `Parameters_*.m` files and wire them into `ScenarioConfig` for full scenario parity.

Generated on: 2025-10-26 (refreshed)
