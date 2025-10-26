## Paper-model MATLAB code — differences vs LLNL baseline (file-by-file)

This document highlights differences between the derivative `paper_model` code and the original LLNL baseline (`LLNL_model`). It lists per-file changes, line-level notable edits (concise, not full dumps), improvements, regressions/risk areas, and recommended next actions.

Scope and method
- I compared the files you provided from `paper_model/` against the previously inspected `LLNL_model/` baseline and recorded functional and stylistic differences. Where the paper variant provides alternate implementations (pump-driven vs pressure-driven), both are discussed.

Top-level summary (high level)
- The `paper_model` branch is an evolved, research-oriented derivative of the LLNL baseline with explicit support for pump-driven transfers, more parameter sets for different scenarios, enhanced plotting and save options, and utilities for logging and results export.
- Key additions: pump control logic (`LH2Control_pump.m`), pump-oriented simulate (`LH2Simulate_Pump.m`), many scenario parameter files (Parameters_*.m), improved plotting (`plotLH2Data.m` now supports saving & nicer layout), and helpers for result saving and logging (`SaveResultsFunction.m`, `UpdateXLSLog.m`). README updated for the paper. Main control script `MAIN.m` provides scenario selection and orchestration.
- Persisting technical debt: heavy REFPROP dependence, extensive use of `evalin('base',...)` and global variables, duplicated helper functions across files, and some occurrences of `clear all` / `close all` in functions which can make reuse and unit testing harder.

Files added/changed and important differences (file-by-file)

- `MAIN.m`
  - Replaces `runNominal.m` as entry point for the paper experiments.
  - Adds a `Case` selector to choose different scenarios (pressure-driven, pump-driven, Topfill/Bottomfill, etc.). This is an improvement: single script to run multiple experimental configurations.
  - Adds numerous orchestration options (SavePlots, SaveResults, WriteXlsTxt, UpdateLog, PCShutdown) and integrates `CreateXLSTXT`, `UpdateXLSLog`, `SaveResultsFunction` and `CreateXLSTXT` hooks.
  - Improvement: more reproducible experiment management and integration with logging/plot saving.

- `Parameters_*.m` (many files: `Parameters_MainToOnboard_PressDiff.m`, `Parameters_MainToOnboard_Pump.m`, `Parameters_Original.m`, `Parameters_TrailerToMain_*` variants)
  - These are multiple alternate parameter sets tailored to specific scenarios (press-diff vs pump, different tank sizes and geometries). They generalize the original single `inputs_TrailerToDewar.m` into scenario-specific configs. This is an important architectural improvement: separating scenario parameters from solver logic.
  - Line-level: many numeric constants differ (volumes, radii, heat-transfer coefficients, vent sizing, pressure setpoints). Example: `LH2Model.VTotal1` and `VTotal2` are different in several files to represent different tank sizes for the paper experiments.

- `LH2Control.m` (paper_model version)
  - Signature changed: now accepts `ETTVentState` as input and returns it in `U.ETVentState` (explicit vent-state feedback). The original LLNL `LH2Control` used `evalin('base','LH2Model')` and global state; paper code still uses `evalin` but passes vent state explicitly: improvement toward clearer control IO.
  - Fill regime breakpoints and multipliers changed. Examples:
    - In LLNL baseline, slow/fast thresholds were `hL2 < 0.05*P.H`, `hL2 < 0.30*P.H`, etc. In paper `LH2Control.m`, fast and reduced-fast and topping ranges are wider (e.g., 0.05 -> 0.70 -> 0.85) reflecting different fill strategy tuning.
    - `lambdaE` (fill valve fraction) values changed (e.g., 0.5/1/1/0.7 in paper vs earlier 0.5/1/1/(1-ET_fill_complete) style in baseline). The paper version consistently multiplies by `(1-ET_fill_complete)` to stop filling once ET full — clearer behavior.
  - Hysteresis widened: getSTVentState uses ±10% thresholds in paper vs ±5% in baseline in some variants — reduces chatter. This is an explicit robustness improvement.

- `LH2Control_pump.m` (new file)
  - New pump-specific control logic. Signature adds `ETTVentState` and uses pump-related parameters `P.PumpMassTransferSlow`, `P.PumpMassTransferFast` to scale `lambdaE`. This file implements pump-driven mode control and is the main functional addition for pump scenarios.
  - Improvement: clean separation of pump-vs-pressure control logic.

- `LH2Simulate.m` (paper_model)
  - Structure largely mirrors LLNL `LH2Simulate`, but with many code paths updated for scenario management and with preallocations and some changed REFPROP call patterns (in some places `refpropm('U','T',P.Tv10,'D',P.rhov10,'PARAHYD')` vs earlier `Q=1` calls). That can be a minor behavior change in initial internal energy setup.
  - Additions: extra global flags added (`ET_vent_complete`, `Process_complete`), and the script starts with `close all; clear all; clc;` — stronger cleanup but breaks interactive workflows and can slow repeated calls. The `waitbar` uses `waitbartime` by default (paper adds a custom waitbar wrapper).
  - In the provided `paper_model` summary, many internal lines were omitted; the visible differences show improved scenario gating and explicit handling for pump vs pressure modes (the `MAIN.m` chooses `LH2Simulate_Pump` for pump cases). The structure preserves event-driven vent switching.
  - Risk: `clear all` inside a function removes breakpoints/compiled state and can cause unwanted side effects; consider `clearvars` local to function instead.

- `LH2Simulate_Pump.m` (new/modified)
  - A pump-oriented variant of the simulation with adjustments to the initial state and control logic. It contains the same physics core but adjusted initializations and control variable handling for pump-driven transfers and additional process flags (`Process_complete`, `ET_vent_complete`, etc.).
  - Improvement: dedicated solver path for pump-driven transfers so experiments can run without mixing control styles.

- `gasFlow.m` (paper_model)
  - Function is unchanged in algorithm (choked/nonchoked logic preserved). Line-level: identical formulas; simply duplicated in `paper_model` for project separation.
  - Neutral: duplication still exists but consistent.

- `vaporpressure.m` (paper_model)
  - Two versions appear in repository historically: an earlier `Tv`/`rhov`-based polynomial and REFPROP branch, and a later `uv`+`rhov`-based function (internal-energy + density). The paper `vaporpressure.m` contains the `uv,rhov` version (robust to going from internal energy output of the solver). This matches `LH2Simulate` which uses internal energies as state variables — alignment is an improvement.
  - Line-level: the paper function attempts `quality=refpropm('q','D',rhov,'U',uv,'PARAHYD')` and truncates `uv` and `rhov` on REFPROP convergence failure (same pattern as baseline, but error messages and truncation differences visible). This is an engineering tradeoff (practical but ad-hoc).

- `plotLH2Data.m` (paper_model)
  - Signature is `plotLH2Data(data,save,path)` and the function supports saving figures and uses nicer default plotting styles (line width, colors, normalized figure size). It also prints summary metrics like EffectiveTransfMass, UsedMass, VentingT1/T2, RelativeVenting at the top — a useful addition for rapid analysis.
  - Improvement: better plotting UX and save-to-path support.

- `SaveResultsFunction.m` and `UpdateXLSLog.m`
  - New utilities for saving `.mat` files and writing summary metrics to Excel logs. `UpdateXLSLog` implements scenario/method classification and writes a results log spreadsheet. These are paper-specific experiment management utilities — clear improvement for reproducibility.

- `CreateXLSTXT.m`, `CreateXLSTXT` referenced in `MAIN.m` (not attached here)
  - The paper repo includes utilities to export main results into text/XLS formats; useful for reproducible figures and supplementary materials.

- `waitbartime.m`
  - Present in `paper_model` as a custom waitbar providing a better UX and optional timing. This is a usability improvement.

- `Data_extraction.m` (paper_model)
  - The paper folder contains its own `Data_extraction.m` (attachment empty in the dataset you provided). In the LLNL baseline, `Data_extraction` computes derived variables and uses a `vaporpressure_bis` helper (which was missing). The paper variant likely adapts extraction for paper-specific outputs: `EffectiveTransfMass`, `UsedMass`, `Boiloff_ST/ET` and so forth (these fields are referenced by `plotLH2Data.m`). If the file is empty or missing content, that is a regression: ensure `Data_extraction.m` is present and consistent with `plotLH2Data` expectations.

Line-level notable differences and examples

- Control and hysteresis
  - Baseline `LH2Control` used ±5% thresholds for vent gating in `getSTVentState`; paper `LH2Control.m` widens this to ±10% in some variants. Lines changed: the conditional threshold arithmetic around `threshold+0.05*threshold` -> `threshold+0.1*threshold`.
  - Example (paper): getSTVentState uses `if p1 < threshold+0.1*threshold` (wider deadband) vs baseline `if p1 < threshold+0.05*threshold`.

- Fill fractions and lambdaE scaling
  - Paper code consistently scales `lambdaE` by `(1-ET_fill_complete)`, explicitly stopping filling when ET full. LLNL code sometimes used `U.lambdaE = (1-ET_fill_complete)` only in specific regimes; paper code applies it more uniformly. This is a safer control policy.

- Pump integration
  - Pump path introduces `P.PumpMassTransferSlow` and `P.PumpMassTransferFast` parameters and uses their ratio to scale `lambdaE` in `LH2Control_pump.m`. This appears in lines like `U.lambdaE = P.PumpMassTransferSlow/P.PumpMassTransferFast*(1-ET_fill_complete)`. This is a new control input and an explicit improvement enabling pump-mode tuning.

- File initialization
  - Paper `LH2Simulate` and `LH2Simulate_Pump` begin with `close all; clear all; clc;` while LLNL baseline did not `clear all` in all versions. This is a functional difference: paper attempts to avoid leftover state but at the cost of clearing workspace and possibly breaking caller state.

- Diagnostic outputs
  - Paper `plotLH2Data` exposes scalar summary metrics (TransferredMass, UsedMass, VentingT1/T2, RelativeVenting) via top-of-file prints and uses percent axes for tank fill. These are new, line-level additions to support paper figures and metrics calculation.

Improvements (what the paper code adds or does better)

1. Pump mode support: clear separation of pump vs pressure-driven control (`LH2Control_pump`, `LH2Simulate_Pump`) and pump-specific parameters make the repository suitable to compare modes and to reproduce paper results.
2. Scenario parameterization: many `Parameters_*.m` files allow easy switching of experimental conditions without editing solver core.
3. Plotting & reporting: `plotLH2Data` improved visuals, saving support, and key scalar metrics printed for quick assessment or manuscript figures.
4. Logging & exports: `SaveResultsFunction` and `UpdateXLSLog` enable experiment provenance and spreadsheets of results for analysis.
5. Robustness: some hysteresis choices widened (±10%) to reduce valve chatter, and the control code more consistently respects `ET_fill_complete` to avoid overfilling.

Regressions / risk areas

1. `clear all` / `close all` usage inside simulation functions — useful for fresh environment but problematic for automation, unit tests, and refactoring into Python modules. Replace with `clearvars -except` or avoid clearing at function scope.
2. REFPROP dependence remains a hard requirement — no CoolProp fallback provided here. Many refpropm calls rely on the REFPROP DLL; this limits portability and CI testing unless a mock backend is introduced.
3. Global variables and `evalin('base',...)` calls persist. While `LH2Control` now accepts vent-state input, many files still use base-workspace coupling. This makes parallel runs or programmatic calls fragile.
4. Duplication persists — `gasFlow`, `cylVToH`, and `vaporpressure` variants appear in multiple places (standalone and nested). Risk of inconsistent updates.
5. Missing / inconsistent helpers: `Data_extraction.m` in `paper_model` appeared empty in the provided snapshot; if truly missing, plots and logging will fail. Similarly, `vaporpressure_bis` used in baseline `Data_extraction` needed reconciliation; check which helper the paper extraction expects.

Recommendations (short-term)

1. Consolidate property access behind an adapter (Paper plan's CoolProp/REFPROP abstraction). Add a mock backend to run tests without REFPROP.
2. Remove `clear all` from `LH2Simulate*` and prefer localised clear or no-clear to improve iterative development and unit testing.
3. Consolidate duplicate helper functions into a single utilities module and reference it from both simulate and post-processing.
4. Confirm and (if needed) re-create `Data_extraction.m` in `paper_model` to match the outputs expected by `plotLH2Data` (fields like `EffectiveTransfMass`, `UsedMass`, `Boiloff_ST`, `Boiloff_ET`).
5. Add a very small regression test (MATLAB script) that runs one short scenario with a mocked property backend to verify integrator/IO plumbing after refactors.

If you want, I can now:
- (A) Produce a one-to-one mapping of LLNL core functions -> paper_model files to create a porting checklist for the Python plan (recommended next step), or
- (B) Begin consolidation work by extracting the common helper functions (`gasFlow`, `cylVToH`, `vaporpressure`) into a single `utils` file in `paper_model` and update call-sites. (I can implement and run basic static checks.)

Tell me which follow-up you prefer and I will proceed.
