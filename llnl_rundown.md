## LLNL MATLAB LH2 Transfer Code — File-by-file rundown

This document summarizes every file in the `LLNL_model` directory, lists each file's purpose, key attributes, important implementation notes, and any issues or dependencies to watch for. At the end is a short overall summary and recommended next steps.

---

### `runNominal.m`
- Purpose: top-level script to run the full simulation workflow.
- Key actions: calls `inputs_TrailerToDewar` to set parameters, runs `LH2Simulate` to run the model, calls `Data_extraction` to post-process and then `plotLH2Data` to visualize. Writes a large `output_bottomfill.txt` with many fields sampled every 10 time steps.
- Notes:
  - Good single-command workflow for reproducing nominal runs in MATLAB.
  - Relies on global workspace behavior (functions call `evalin('base',...)` and expect `LH2Model` in base workspace).

### `inputs_TrailerToDewar.m`
- Purpose: initialize model parameters and initial conditions in the `LH2Model` struct placed in base workspace.
- Key attributes:
  - Geometry: trailer (ST) horizontal cylinder and Dewar (ET) vertical cylinder. Volumes, radii, areas, lengths set.
  - Thermophysical constants and correlations: critical T, p, liquid density, specific heats, thermal conductivities, viscosities, vapor gas constants, gamma, etc.
  - Grid sizing: `nL1`, `nV1`, `nL2`, `nV2` etc. used by solver discretization in `LH2Simulate`.
  - Vent, valve, pipe, transmission-line, vaporizer and solver-related parameters (e.g., `tFinal`, `relTol`).
- Dependencies: calls `refpropm` (REFPROP) to get initial densities and uses polynomial fits to REFPROP outputs.
- Notes:
  - Many magic numbers and polynomial fits are embedded; pay attention to units (some REFPROP calls use g/L etc.).
  - Comments caution about solver stiffness and recommended tuning (e.g., reduce valve areas or increase time constants).

### `LH2Simulate.m`
- Purpose: core ODE simulation implementation. Integrates mass/energy balances, handles events and mode switching, returns `data` structure with time series and many diagnostic fields.
- Key structure:
  - Builds initial state vector `x0` containing liquid masses, internal energies per grid cell, vapor masses/energies, surface temperature, transfer flows, wall temperature, and many diagnostic slots.
  - Defines nested function `LH2dxdt(P,t,x)` which computes derivatives (mass and energy balances, heat transfer, condensation, pdV work, valve flows, vaporizer dynamics, etc.). Uses REFPROP (`refpropm`) extensively for densities, temperatures, enthalpies, transport properties and qualities.
  - Implements event-driven stepping: uses `ode15s` with an `Events` function `VentEvents` to stop when ET vent conditions change; then toggles vent state and continues (discrete-event handling via re-starting ODE solve windows).
  - Contains local helper functions: `dsqrt` (directed sqrt), local `gasFlow` (choked/non-choked flow), `cylVToH` (cylinder volume-to-height), `vaporpressure` (internal-energy+density-based vapor pressure), and `VentEvents` which defines event roots for vent thresholds.
- Key physics / numerics notes:
  - Uses polynomial correlations (from REFPROP) for liquid density vs internal energy and enthalpy of vaporization vs film temperature.
  - Heat transfer uses laminated boundary-layer discretization (grids `l12_*`, `l_*`) for vapor and liquid layers.
  - Handles two-phase vs supercritical behavior: uses `refpropm` quality checks and may fall back to truncated inputs when REFPROP doesn't converge.
  - `ode15s` used with event detection; integration is restarted when vent state toggles to enforce discrete valve logic.
  - Many diagnostic variables duplicated and appended to state vector so they are stored in output arrays.
- Implementation issues/observations:
  - Heavy reliance on REFPROP MEX / DLL (`refpropm`) which must be installed and on path; otherwise many calls will fail.
  - There is internal duplication: `gasFlow` and `cylVToH` functions also exist as separate files. The nested definitions inside `LH2Simulate` shadow external copies when `LH2Simulate` runs.
  - Frequent use of `evalin('base',...)` / `assignin('base',...)` and global flags (`ETTVentState`, `ET_fill_complete`, etc.). This makes pure-function testing harder.
  - Long function with many responsibilities; modularization would reduce coupling and improve testability.

### `LH2Control.m`
- Purpose: simple state machine / controller that sets valve operating fractions depending on ET liquid height and ST/ET pressures.
- Key behavior:
  - Decides `lambdaE` (fill valve opening fraction), `lambdaV` (vaporizer valve), and `STVentState` based on `hL2` and pressure setpoints.
  - Embeds small helper functions: `getSTVentState`, `getETVentState`, `getVaporizerValveState` implementing hysteresis-like behavior.
- Notes:
  - Uses `evalin('base','LH2Model')` to access parameters.
  - Control logic is simple and rule-based (thresholds ±2–5%).

### `Data_extraction.m`
- Purpose: post-process simulation output (often the `nominal` struct) to compute derived quantities and prepare arrays for plotting and saving.
- Key operations:
  - Computes liquid volumes from masses and polynomial density correlation (same polynomial used elsewhere) and ullage volumes and vapor densities.
  - Calls `vaporpressure_bis(...)` in this file (note: this symbol is referenced but not present in the repository — see "Issues" below).
  - Calls `cylVToH` (the separate file copy) to convert horizontal cylinder liquid volume to height.
  - Computes boiloff accumulation `Boiloff_ET` as time integral of vent flows.
- Notes / issues:
  - `Data_extraction.m` uses `vaporpressure_bis(...)` but the repository contains `vaporpressure.m` and the nested version inside `LH2Simulate`. There may be a missing `vaporpressure_bis` helper or alternate version; investigate to ensure consistent usage.

### `plotLH2Data.m`
- Purpose: create a set of diagnostic plots using the `data` struct returned by `LH2Simulate`.
- Key plots produced:
  - Liquid levels and transfer line mass flow.
  - Temperatures for liquid, surface, vapor, and wall in ST and ET.
  - Pressures and internal energies.
  - Masses and densities.
  - Venting fluxes and heat transfer breakdowns (several grouped subplots and diagnostic legends).
- Notes:
  - Expects `LH2Model` to be available in base workspace (uses `evalin`).
  - Legend and axis choices reflect the original LLNL plotting conventions.

### `cylVToH.m`
- Purpose: for a horizontal cylinder (ST/trailer) convert a liquid volume into an equivalent fluid height (geometry math).
- Implementation details:
  - Computes cross-section area and reduces the problem to solving a nonlinear equation for the half-height parameter `x`.
  - Uses a Newton iteration with finite expressions for derivative `supd` until relative error < 1e-4.
  - Returns `H = R ± x` depending whether the filled fraction `s` is above half the cross-sectional area.
- Notes:
  - A nearly identical `cylVToH` function is also defined as a nested function inside `LH2Simulate` (duplication).

### `gasFlow.m`
- Purpose: compute choked or non-choked mass flow rate through an orifice/orifice-like geometry between two pressures P1 and P2 for a perfect-gas exponent `gamma` and density `rho`.
- Implementation:
  - If P1 < P2, it returns negative of swapped arguments (so flow sign reverses correctly).
  - Determines choked condition using threshold = ((gamma+1)/2)^(gamma/(gamma-1)).
  - Choked formula: mdot = CA * sqrt(gamma * rho * P1 * (2/(gamma+1))^((gamma+1)/(gamma-1))).
  - Non-choked formula implemented using isentropic relationships.
- Notes:
  - There is a second copy of `gasFlow` inside `LH2Simulate` (duplication).

### `vaporpressure.m`
- Purpose: compute vapor pressure (Pa) for a given film temperature and vapor density; handles two-phase and supercritical regimes.
- Implementation details:
  - Uses a polynomial `T_sat(rhov)` to get a saturation temperature curve from density (polynomial given in code).
  - If provided film temperature `Tv` > `T_sat` (supercritical branch), calls REFPROP `refpropm('P','T',Tv,'D',rhov,'PARAHYD')*1e3` to compute pressure.
  - If two-phase branch, uses a polynomial cubic in `Tv` (coefs are embedded) to compute `pv` (returned in Pa).
- Notes:
  - `LH2Simulate` contains a different vaporpressure implementation (which uses quality tests and `refpropm` paths) — there are multiple variants in the codebase. Careful consolidation is recommended.

### `README.md`
- Purpose: repository-level documentation describing LH2TS, references, files list, and REFPROP dependency.
- Key notes:
  - Documents that the code was developed/tested with MATLAB R2013b and REFPROP 9.1.
  - Advises on REFPROP linking and DLL placement for MATLAB.

### `LICENSE.md`
- Purpose: the repository includes the original license text. The code is under the NASA Open Source Agreement (NOSA) v1.3 (see file for terms and obligations).
- Important: NOSA has obligations for redistribution and must be respected when deriving new works; check compatibility with any chosen Python license if this code is used as a behavioral reference.

---

## Cross-file notes and issues observed

- REFPROP dependency: The code extensively uses `refpropm` (the REFPROP MATLAB wrapper) for EOS and transport properties. REFPROP must be correctly installed and visible to MATLAB for the simulations to run. Several catch/fallbacks exist to mitigate non-convergence but they are ad-hoc (truncation of uv/rhov values).
- Duplication: `LH2Simulate` contains local copies of several helper functions (`gasFlow`, `cylVToH`, `vaporpressure`) which are also present as standalone `.m` files. This duplication can cause maintenance problems and subtle inconsistencies (e.g., `Data_extraction` calls `vaporpressure_bis` which is not part of this folder).
- `vaporpressure_bis` missing: `Data_extraction.m` uses `vaporpressure_bis(...)` while the repo contains `vaporpressure.m` and the nested version inside `LH2Simulate`. Confirm whether `vaporpressure_bis` is a separate helper in another branch or was intended to alias `vaporpressure`.
- Global state / base workspace reliance: controller and parameter passing relies on `evalin('base',...)` and `assignin('base',...)` plus global flags (`ETTVentState`, etc.). These practices make reuse (for example, calling functions from a unit-test harness or a Python port) harder and increase risk for hidden coupling.
- Event handling style: the ODE integration restarts on vent events (`VentEvents` returns roots and integration is restarted with toggled vent state). This is a valid approach but requires careful event definitions; the code uses `Events` and `isterminal` settings for robust switching.

## Overall codebase summary

This LLNL MATLAB codebase implements a detailed lumped-parameter, physics-informed model of liquid hydrogen transfer between a horizontal trailer (ST) and a vertical Dewar (ET). It simulates liquid and vapor control volumes in each tank, uses REFPROP for thermophysical properties, and includes heat transfer, condensation/evaporation, valve/vent flow (including choked flow), a vaporizer tube model, and event-driven control (venting thresholds, fill-stage controller). The `LH2Simulate` routine is the core numerical driver (uses `ode15s`), and `runNominal.m` ties the workflow together. The repository includes plotting and postprocessing utilities.

Strengths:
- Physics-rich model with many practical details (transport properties, condensation physics, film temperature correlations, vent logic).
- Event-driven integration and diagnostic outputs collected for detailed post-analysis.

Weaknesses / modernization opportunities:
- Heavy REFPROP dependence and many inline polynomial fits which should be centralized for maintainability.
- Duplication of helper functions; consolidate single-source helpers.
- Global workspace reliance and nested long functions — refactor to smaller pure functions and explicit input/output for easier testing and for porting to Python.
- Missing helper reference `vaporpressure_bis` referenced from `Data_extraction.m` — needs reconciliation.

Recommended next steps if you plan to port or reuse this code as the basis for a Python project:
1. Consolidate thermophysical property calls behind a single adapter layer (e.g., a `Properties` module) and add a mock backend for tests. That will make switching REFPROP/CoolProp straightforward.
2. Remove base-workspace globals where possible. Make `LH2Simulate` accept a parameter struct and return outputs without writing to base workspace.
3. Consolidate duplicated helpers (`gasFlow`, `cylVToH`, `vaporpressure`) into single utilities and add unit tests for each.
4. Locate or recover `vaporpressure_bis` or update `Data_extraction.m` to call the canonical `vaporpressure` implementation.
5. Add small regression scenarios and unit tests for core numeric building blocks (cylinder geometry, choked-flow formula, vapor-pressure approx.) before major refactor.

---

If you want, I can now:
- (A) Create a small checklist / mapping to port each function to Python modules (properties, geometry, unitops, control, integrator, IO, viz) matching your project plan; or
- (B) Start extracting the minimal, well-contained pieces (e.g., `gasFlow`, `cylVToH`, `vaporpressure`) into separate testable files or Python equivalents.

Tell me which follow-up you prefer and I'll proceed.
