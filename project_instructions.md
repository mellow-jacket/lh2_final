Python LH₂ Transfer Modeling Framework — Implementation Plan

Project codename: h2_delivery_model

0) Executive summary

We will build a modern, open, and extensible Python framework that replicates and then extends the MATLAB-based LH₂ transfer models (LLNL’s LH2Transfer and the UCI/CEI derivative work) with modular unit operations, configurable geometries, event/mode logic, and traceable artifacts. The work program below converts the “challenge surface” you outlined into an actionable, staged plan that specifies what to build (APIs, modules, tests, artifacts) and how we will know it’s done (acceptance evidence and reproducibility), while avoiding prescriptive solver formulas or parameter values.

Alignment with the reference paper and repos. The CEI paper adds four capabilities over the LLNL baseline—pump model, line/valve pressure losses, supply‑tank wall thermal mass & vent, and geometry options—and formalizes metrics (relative venting, compression ratio, SOC) we will use for acceptance. The figures and equations in §2–§3 of the paper define the qualitative patterns our Python results must reproduce (e.g., vent clustering near high SOC; non‑monotonic maximum transferred mass vs flow rate) without enforcing numeric identity. We will wire those figures into our validation suite as traceable scenarios. 

h2_delivery_model

1) Purpose, scope, and “done” (acceptance-only)

Purpose. Deliver a Python package that can simulate pressure‑driven and pump‑driven LH₂ transfers between two tanks with event logic (vent/MWP/SOC), pluggable thermophysical properties, configurable geometries, and provenance‑rich artifacts (CSV/Parquet and figures) suitable for reproduction of the qualitative results and event sequences in the paper (e.g., Fig. 1 schematics on p.4; venting surfaces in Figs. 7–8; SOC effects in Figs. 9–10; pump vs pressure‑driven comparison in Fig. 11; MTM behavior in Fig. 12). 

h2_delivery_model

Non‑goals in v1. No dormancy modeling, no para↔ortho conversion, no economics/optimization, no prescription of solver formulas or parameter values. (The paper also assumes 100% parahydrogen and neglects spin conversion over short transfer times.) 

h2_delivery_model

Success criteria (evidence we will collect).

Qualitative parity: Across curated scenarios, replicate the shapes and trends reported (e.g., venting rises with flow rate and initial receiving‑tank pressure; venting clusters near high SOC; MTM vs flow is non‑monotonic). Tie each figure to a versioned scenario file. 

h2_delivery_model

Event correctness: Mode switches (vent open/close, MWP trips, pump enable/disable, choked↔unchoked venting) occur in a causal order without conservation breaks across state resets.

Reproducibility: Deterministic trajectories/artifacts (within declared tolerances) per OS/CPU when configuration & property backend are fixed.

Traceability: Every artifact contains scenario hash, property backend/version, geometry IDs, and code version; figures regenerate from scenarios.

Envelope adherence: Clear warnings or halts when property calls are out of declared validity ranges.

2) Reference alignment & IP guardrails

LLNL baseline (MATLAB): Code exists as LH2Transfer under the NASA Open Source Agreement 1.3; parts of that code (and embedded assumptions) guide behavior, terminology, and data structures we mirror—but we will not copy code. We will attribute and keep license boundaries clear. 
GitHub
+1

Derivative repo (Albert‑Gil / MATLAB): GPL‑3.0; used as a behavioral reference and for figure replication targets; again, no code reuse. 
GitHub

License notes: NASA OSA 1.3 is OSI‑approved but has terms that can complicate code mixing; we will avoid incorporating NOSA‑licensed code and implement a clean‑room Python design with permissive licensing (e.g., BSD‑3/MIT). 
Open Source Initiative
+1

3) Thermophysical properties & units (plug‑ins, not prescriptions)

Backends (pluggable):

CoolProp (default, open‑source MIT; supports ParaHydrogen): easy to distribute in CI; our Para‑H₂ baseline will target CoolProp’s fluid model for validation runs. 
GitHub
+1

REFPROP (optional, licensed): higher‑fidelity real‑gas data; Python wrapper ctREFPROP available; we will support it behind a feature flag and mark CI tests that require it as “optional.” Distribution of public online apps requires NIST permission per license. 
GitHub
+2
PyPI
+2

Units: Adopt Pint for a project‑wide unit registry and IO validation (bar vs Pa, W vs kW, etc.). Configs declare units explicitly and are validated on load. 
Pint

Property envelope & uncertainty: Each backend declares a valid state envelope; we log every out‑of‑envelope request and surface it as a test failure. Backends can report a (documented) version string captured in artifacts.

4) System architecture (module boundaries & APIs)

We architect for clarity, extensibility, and traceability, not for a single solver choice.

Package layout (top‑level namespaces):

h2_delivery_model.properties — backend adapters (CoolProp/REFPROP), with an abstract PropertyBackend that exposes: (T,P)->ρ,h,cp,μ,k,γ,Z,q_vap, phase queries, backend/version, and envelope checks.

h2_delivery_model.geometry — shape models (vertical cylinder, horizontal cylinder, sphere) to compute level↔volume, wetted area, interface area and head; exactly the variants described in §2.1.4. 

h2_delivery_model

h2_delivery_model.unitops — composable unit operations: Tank, Pump, Valve, Pipe, Junction. Tanks expose ports and internal control volumes (liquid, vapor, interface, wall lump(s)).

h2_delivery_model.events — declarative EventRule with guards (e.g., P>=MWP, SOC>=max, choked threshold) and effects (mass/energy removal via vent, pump enable/disable, mode flag flips) with hysteresis fields to prevent chatter.

h2_delivery_model.integrators — thin adapters around numerical engines exposing a unified step/event API. We will support a pure‑Python ODE route and a DAE route (see §6), but which solver is used remains a configuration choice.

h2_delivery_model.sim — the orchestrator that wires flowsheets, state vectors, event loop, and state re‑projection at discrete jumps.

h2_delivery_model.scenario — Pydantic‑validated YAML/JSON configs with units, geometry, fittings, initial states, and an operating schedule (e.g., pump rate profile).

h2_delivery_model.io — readers/writers to Parquet/CSV and attribute‑rich NetCDF/xarray datasets; provenance stamping. 
xarray

h2_delivery_model.viz — plotting helpers to regenerate figures comparable to the paper’s conventions (axes, units) in headless CI.

h2_delivery_model.testing — invariant checks, scenario regression harness, property mocks.

Data model & artifacts.

Time series: stored as xarray.Dataset with labeled dims: time, tank ∈ {1,2}, phase ∈ {L,V}; dataset attrs store scenario hash, code version, property backend & version. 
xarray

Events: structured log (JSON Lines / CSV) containing event name, time, guard values, and mass/energy deltas; persisted with provenance.

Parquet: write per‑run tables with key‑value metadata for external tools; store pointers to scenario & code version. 
mungingdata.com

API sketch (pseudocode):

cfg = Scenario.from_yaml("pump_transfer_sweep.yaml")   # units validated (Pint)
sim = Simulator(cfg, properties=CoolPropBackend("ParaHydrogen"))
run = sim.run()                                        # time series + events + KPIs + figures
run.save("artifacts/2025-05-01T12-30Z")

5) Physics scope and extensibility (no formulas prescribed)

We mirror the conceptual partitions in the reference work (§2): three control volumes per tank (liquid, vapor, mass‑less interface) plus wall thermal mass; eventful transfers driven by pressure difference or by a pump; line losses for valves and pipes; and vent flows that transition between choked and non‑choked based on real‑gas compressibility—exactly the pieces emphasized in the methodology and Fig. 1 schematics (p.4). 

h2_delivery_model

Key extension seams (v1 interfaces in parentheses):

Pump injection (unit op): enthalpy rise tied to work/efficiency; decouples flow from ΔP as described in §2.1.1. 

h2_delivery_model

Line/valve losses (unit ops): K and Darcy–Weisbach abstractions (§2.1.2). 

h2_delivery_model

Supply‑tank wall mass & vent (tank extensions) (§2.1.3). 

h2_delivery_model

Geometry variants (geometry module) (§2.1.4). 

h2_delivery_model

Interface dynamics (model kernel): expose the condensation‑blocking mode switches highlighted in §2.2.3 without encoding specific closures here. 

h2_delivery_model

6) Numerical integration & event handling (choices, not mandates)

We design a pluggable integration layer with two families:

ODE path (SciPy): scipy.integrate.solve_ivp with event functions for vent/MWP/SOC/choking roots, acknowledging root‑finding semantics (sign change per step; multiple zero crossings may be missed unless step control is tuned). 
SciPy Documentation

DAE path (SUNDIALS via scikits‑odes / Assimulo / scikit‑SUNDAE): expose IDA/CVODE integration for stiff systems and index‑1 DAEs, with built‑in event support and robust step rejection around mode switches. We will keep the Python API stable while allowing users to switch backends. 
PyPI
+2
GitHub
+2

Discrete events & re‑projection. Our EventEngine applies conservative state updates (e.g., vent removes m & energy per vent model; pump toggles setpoints), then re‑projects invariants (non‑negative masses, sum constraints). A hysteresis/deadband field on each rule reduces chatter and Zeno‑like toggling. All event effects are logged for audit.

Determinism: Document thread/process guidelines; suggest process‑based parallelism for sweeps (e.g., joblib loky) to avoid global state in property DLLs. Pin BLAS threads in CI. 
joblib.readthedocs.io

7) Configuration & provenance

Config schema (Pydantic):

Tanks: geometry (shape + dimensions), wall mass options, MWP/SOC limits, vent parameters, initial T/P/SOC.

Network: pipes (L, D, roughness), valves (Cv or K), pump (efficiency envelope or neutral map placeholder), elevations.

Operation: mode selection (pressure vs pump), pump flow profile, vent bands, stopping criteria (e.g., SOC, time).

Backend & reproducibility: property backend/version pin; random seed (if used); solver backends & tolerances (names only; no values prescribed).

Every run generates a manifest stamped into artifacts: Git commit, scenario hash, backend version, OS/CPU.

8) Validation strategy (what we’ll reproduce and how we’ll show it)

We curate reference scenarios that recreate the narrative experiments:

Pressure vs pump differences in pressure trajectories, controlled vs uncontrolled flow, and temperature trends—compare to Figs. 2–3 (§3.2). 

h2_delivery_model

Pump‑driven dynamics and vent bursts around high SOC—compare to Fig. 4 (§3.3.1). 

h2_delivery_model

Flow‑rate extremes (“slow” vs “fast”) and heat‑flow sign changes—compare to Figs. 5–6 (§§3.3.2–3.3.3). 

h2_delivery_model

Relative venting surfaces vs flow rate and initial pressure (Fig. 7) and vs inlet LH₂ temperature (Fig. 8). 

h2_delivery_model

SOC dependence with variable and fixed transferred mass (Figs. 9–10). 

h2_delivery_model

Pump vs pressure‑driven venting comparison (Fig. 11). 

h2_delivery_model

MTM non‑monotonicity vs flow rate (Fig. 12) and the coldest final LH₂ at intermediate flows—capture the shape and turning points qualitatively. 

h2_delivery_model

Mass/energy closure: The paper reports 0.9–1.7% balance closure for simulated cases; our acceptance will require comparable closure envelopes (documented in test plans) while not prescribing solver tolerances. 

h2_delivery_model

Metric definitions: Implement Relative venting, Compression ratio, and SOC exactly as defined in §2.3 (Eqs. 35–36) to ensure apples‑to‑apples comparisons in KPIs and plots. 

h2_delivery_model

9) Software engineering: quality, packaging, CI/CD

Modern packaging: single‑source pyproject.toml metadata per PEP 621; build via Setuptools or Hatch; wheels for Linux/macOS/Windows; optional extras refprop, dev, docs. 
Python Packaging
+2
setuptools.pypa.io
+2

Code quality & hooks: adopt Ruff (linter/formatter), pytest, mypy (opt‑in), and pre‑commit hooks (run ruff/pytest‑quick). 
Astral Docs

CI: GitHub Actions matrix across OSes and Python versions; headless Matplotlib backend; two test tiers:

Core: CoolProp‑backed scenarios + property mocks (always on).

Optional: REFPROP‑backed scenarios (skipped in public CI unless REFPROP is available). REFPROP tests document local setup (ctREFPROP) and DLL pathing quirks. 
GitHub

Logging/telemetry: structured logs via structlog with JSON output option; event audit trail shipped as artifacts. 
structlog

Artifacts: Parquet/CSV + PNG/PDF; each embeds metadata (scenario hash, backend/version). 
mungingdata.com

10) Performance & scaling

Hotspots: property evaluations, event detection near switching, tight Python loops.
Mitigations: vectorize where possible; small memoization windows for repeated property calls; optional Numba on pure‑Python kernels that don’t call into external DLLs; scenario sweeps via process‑based parallelism (joblib/loky). 
joblib.readthedocs.io

11) Risk register & mitigations (selected)
Risk	Indicator	Mitigation
Backend drift (CoolProp/REFPROP updates shift properties)	Regression mismatches in vent timing/pressure peaks	Pin backend versions in scenarios; include backend version in artifacts; keep dual‑backend tests. 
GitHub

Event instability (vent chatter, simultaneous triggers)	Step rejection storms; negative inventories	Require hysteresis in rules; sequence simultaneous roots deterministically; post‑event re‑projection with invariant checks.
License/CI friction (REFPROP)	Skipped tests; platform‑path errors	Optionalize REFPROP; document ctREFPROP setup; use process isolation; public CI runs CoolProp only. 
GitHub

Divergence from references	Missing vent clustering at high SOC; MTM not visible	Use paper’s scenarios/figures as acceptance templates; cross‑check KPI definitions (§2.3). 

h2_delivery_model

12) Roadmap & milestones (entrance/exit = acceptance only)

M0 — Foundations
Entrance: Repo skeleton, coding standards, unit registry chosen.
Exit: Package imports; config schema validates; unit‑checked logging; CI “core” green.

M1 — Baseline physics (pressure‑driven)
Entrance: Property backend adapter passes envelope tests (Para‑H₂ in CoolProp).
Exit: Two‑tank pressure‑driven scenario runs; event logs show vent bands and pressure trends comparable to paper (Fig. 3 left panel shapes). 

h2_delivery_model

M2 — Pump & line/valve losses
Entrance: Pump and line/valve unit‑ops registered.
Exit: Pump‑driven scenarios match qualitative differences vs pressure‑driven (controlled flow, reduced venting; Fig. 11 trends). 

h2_delivery_model

M3 — Geometry & wall thermal mass parity
Entrance: Geometry module (vertical/horizontal cylinder, sphere); wall mass in both tanks.
Exit: Scenarios demonstrate shape/thermal mass impacts; artifacts stamped with geometry metadata (§2.1.4). 

h2_delivery_model

M4 — Validation portfolio
Entrance: Scenario set for sweeps (flow rate, initial pressure, inlet temperature, SOC).
Exit: Reproduced venting surfaces (Figs. 7–8), SOC effects (Figs. 9–10), MTM curves (Fig. 12) with expected qualitative shapes and clustering near high SOC. 

h2_delivery_model

M5 — Packaging & release
Entrance: Tests green locally; artifacts reproducible.
Exit: Cross‑platform wheels; docs & examples; release with embedded provenance.

13) What we need up front

Decisions:

Default license (recommend BSD‑3/MIT) to keep downstream use open; keep strict separation from NOSA/GPL code. 
Open Source Initiative
+1

Initial property backend (CoolProp ParaHydrogen as default; optional REFPROP). 
coolprop.org
+1

Accepted OS matrix (Windows/macOS/Linux) and Python versions.

Inputs: geometry & valve/pipe data for example scenarios; vent setpoints and MWP/SOC limits for demo cases; agreed plot styles/units.

14) Testing strategy (how we keep it correct)

Unit tests: geometry mappings (level↔volume/area), SOC/MWP logic, vent choked/un‑choked switching, conservation invariants.

Property mocks: “table” backend for deterministic E2E tests; backends exercise envelope violations.

Scenario regression tests: golden artifacts (time series + KPIs) for a small suite of canonical cases; numerical tolerances documented per signal.

Performance tests: simple timing baselines on one or two scenarios to catch regressions.

Cross‑platform smoke tests: short runs in CI on all OSes.

15) Documentation, examples, and education

User guide: four progressive examples that mirror the paper narratives—(i) minimal single‑tank warm‑up, (ii) two‑tank pressure‑driven, (iii) pump‑driven with vent events, (iv) sensitivity sweeps producing surfaces. Each example ships a scenario file + one‑command “regenerate figures” script. 

h2_delivery_model

API reference: autodoc, with short runnable snippets.

Repro kits: requirements.txt/pyproject.toml, backend versions, and a “How we pin backends” section.

16) Interoperability & future bridges (optional, post‑v1)

Optimization frameworks: Provide a data bridge (not tight coupling) to Pyomo/IDAES for steady‑state or coarse dynamic studies later; DAEs are naturally supported there, but our v1 will remain physics‑first. 
pyomo.readthedocs.io
+1

Batch runs: CLI to execute sweep grids and collect KPIs into a Parquet results table with metadata.

17) What the finished product looks like (user view)

CLI (examples):

# Run a pressure-driven transfer reference and regenerate paper-style plots
h2model run scenarios/pressure_driven_reference.yaml --out artifacts/run_001

# Run a sweep (flow rate × initial pressure) and build the relative-venting surface
h2model sweep scenarios/pump_sweep.yaml -o artifacts/sweep_venting


Python API (examples):

from h2_delivery_model import Scenario, Simulator, CoolPropBackend

scn = Scenario.from_yaml("scenarios/pump_driven.yaml")
sim = Simulator(scn, properties=CoolPropBackend("ParaHydrogen"))
run = sim.run()

# KPIs aligned with §2.3 (relative venting, CR, SOC)
print(run.kpis["relative_venting"])
print(run.kpis["compression_ratio"])

# Artifact bundle with provenance (xarray + Parquet + figures)
run.save("artifacts/2025-04-30T18-20Z")


Artifacts:

Time series (xarray/NetCDF) with attrs: scenario_sha, code_version, property_backend, backend_version.

Events log (CSV/JSON) with guard values and effects.

Figures comparable to: venting surfaces (Figs. 7–8), SOC pressure traces (Fig. 10), MTM curve (Fig. 12). 

h2_delivery_model

18) Library choices & best‑practice references (why these)

Properties: CoolProp (open‑source, ParaHydrogen supported; MIT license), optional REFPROP via ctREFPROP with license awareness. 
coolprop.org
+2
GitHub
+2

Solvers: SciPy solve_ivp for ODEs with event roots; SUNDIALS via scikits‑odes/Assimulo/scikit‑SUNDAE for stiff ODE/DAE and robust event handling. 
GitHub
+3
SciPy Documentation
+3
PyPI
+3

Units & data: Pint for dimensional safety; xarray for labeled multidimensional arrays with metadata; Parquet for fast tabular KPIs with custom metadata. 
Pint
+2
xarray
+2

Parallel sweeps: joblib with loky backend (process‑based). 
joblib.readthedocs.io

Packaging & CI: PEP‑621 pyproject; Setuptools; ruff + pre‑commit; GitHub Actions matrix. 
Python Packaging
+2
setuptools.pypa.io
+2

19) Work breakdown (high‑level tasks)

A. Foundations (M0)

Create repo with pyproject.toml (PEP‑621), CI skeleton, ruff/pre‑commit, pytest harness. 
Python Packaging

Add Scenario schema (Pydantic) with units (Pint).

Implement PropertyBackend interface + CoolProp adapter (Para‑H₂). 
coolprop.org

Add logging (structlog) and artifact writer with provenance. 
structlog

B. Baseline physics (M1)

Implement Tank with liquid/vapor/interface/wall lumps; implement pressure‑driven transfer path.

Event engine: MWP/SOC/vent rules with hysteresis; mass‑/energy‑conserving event effects.

First acceptance: pressure‑driven scenario replicates key shape trends vs paper (Fig. 3 left). 

h2_delivery_model

C. Pump & line/valves (M2)

Add Pump, Valve, Pipe unit ops; adopt line loss abstractions.

Add pump‑driven mode with external profile; acceptance versus Fig. 11 trend. 

h2_delivery_model

D. Geometry & wall parity (M3)

Geometry module for vertical/horizontal cylinders and sphere; both tanks get wall thermal mass + vent valves as in §2.1.3/2.1.4. 

h2_delivery_model

E. Validation portfolio (M4)

Scenario sweeps to regenerate venting surfaces and MTM curves (Figs. 7–8, 12). 

h2_delivery_model

KPI calculators (Relative venting, CR, SOC) per §2.3. 

h2_delivery_model

Mass/energy closure checks targeting the paper’s reported closure envelope. 

h2_delivery_model

F. Packaging & Docs (M5)

Sphinx or MkDocs site with examples; versioned release; wheels for 3 OSes.

20) “Edge‑case” playbook (captured in tests)

Empty or near‑empty tanks; high SOC; reverse flow; simultaneous guard triggers.

Vent choking thresholds; transitions with hysteresis; minimum suction pressure for pump; ceiling/floor guards for T/P.

Out‑of‑envelope property calls → warning or halt with explicit diagnostics.

Appendix A — Traceability to the reference paper

Fig. 1 (p.4): mass/energy balance schematics for pressure‑ and pump‑driven cases guide our unit‑op boundaries and ports. 

h2_delivery_model

§2.1.1–2.1.4: define pump, line/valve losses, supply‑tank wall & vent, and geometry options—the feature parity we target in v1. 

h2_delivery_model

§2.3: defines Relative venting, Compression ratio, SOC—our KPI implementations and plot labels use these exact definitions. 

h2_delivery_model

§3 & Figs. 2–12: provide the qualitative acceptance yardsticks we will reproduce (vent clustering near high SOC; pump reduces venting; non‑monotonic MTM vs flow). 

h2_delivery_model

Mass/energy closure (p.7): paper reports 0.9–1.7%—we set comparable closure goals in tests (without prescribing specific tolerances). 

h2_delivery_model

Appendix B — Notes on third‑party libraries and licenses

CoolProp: MIT‑licensed; ParaHydrogen supported; suitable as default in open CI. 
GitHub
+1

REFPROP: accessible via ctREFPROP wrapper; per NIST, public Internet access requires a distribution agreement—so we keep it optional and documented. 
GitHub
+1

SciPy events: event root detection and caveats (missed multiple roots within a single step if not configured). 
SciPy Documentation

SUNDIALS bindings: scikits‑odes (IDA/CVODE), Assimulo, or scikit‑SUNDAE as alternatives if DAE stiffness or event robustness needs increase. 
PyPI
+2
GitHub
+2

Packaging: PEP‑621 pyproject.toml; Setuptools config. 
Python Packaging
+1

Closing

This plan translates the complexities you enumerated into a concrete, staged build with strong traceability to the paper’s design (pump, losses, geometry, wall mass, events) and its acceptance targets (venting surfaces, SOC effects, MTM curves, event timing behavior). It delivers a clean‑room Python codebase with modern packaging, testing, and provenance so researchers can reproduce the paper’s trends and extend the framework safely over time. 

h2_delivery_model

If you’d like, I can draft the initial repository skeleton (pyproject.toml, package layout, config schema, and CoolProp adapter) and the first two validation scenarios next.