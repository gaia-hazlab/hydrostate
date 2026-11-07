# Project instructions

## Purpose and current scope

Develop a physically defensible, reproducible framework coupling DHSVM snow,
vegetation and land-surface hydrology with vadose-zone and MODFLOW groundwater
flow, constrained by observations through data assimilation. Evaluate watershed,
Mount Rainier and western Washington applications, with hourly or finer
hydrologic states and a workstation as the initial compute target.

Hazard applications include earthquake liquefaction, landslides, lateral
spreading, floods, wildfire conditioning and postfire hazards. Hydrologic
conditions must be passed to appropriate hazard models before making claims
about triggering, deformation, inundation or runout.

Hydrostate contains a monthly precipitation/recharge pipeline, a MODFLOW 6
software template, independent synthetic 2D hillslope and 3D valley Richards
experiments, analytical/manufactured verification and sensor CSV operators.
It also provides fixed-grid conservative material updates, complete hydraulic
checkpoints and typed Gaia DataHub interchange for a companion Mechanicalstate
disturbance/recovery runner. These are prescribed synthetic material scenarios,
not mechanics or calibrated damage/healing laws. See
`docs/technical/material_updates.md` and `docs/datahub_contract.md`.
The solvers are not a DHSVM–MODFLOW coupler. The 2D outlet hydrograph is not
mesh-converged; the 3D valley has not established solution mesh convergence.
The linear verification branch does not establish all nonlinear accuracy.
Synthetic sensor comparisons are not field validation or online assimilation.
Read `docs/technical/` and `docs/use_cases.md` for reproducible experiments,
`docs/toy_watershed.md` for measured checks and limitations,
`docs/vadose_design.md` for the target process/coupling contract and
`docs/spatial_parallelism.md` for optimization/decomposition plans.
The project is not yet a calibrated coupled hydrologic or hazard model.
Preserve that distinction in code, documentation and products. Do not
reinterpret smoke-test results or budget closure as scientific validation.

## Physical requirements

- State the governing equations, approximations, boundary conditions and
  applicable regimes for each process. Distinguish an established model
  capability from a proposed implementation or research fork.
- Conserve water across snow/ice, interception, soil, vadose, saturated,
  channel, lake and surface stores included in an application. Include
  precipitation, melt-related changes in ice storage, ET, pumping and all
  external boundary flows in the budget. Internal transfers cancel in the
  coupled budget.
- Report signed budget residuals in water volume and equivalent depth, along
  with a normalized error whose denominator is stated. Set justified acceptance
  tolerances for each experiment; do not hide local errors behind global
  cancellation or accept model termination as evidence of closure.
- Distinguish water content, saturation, matric/pressure head, hydraulic head,
  pore-water pressure, water-table elevation and confined potentiometric head.
  Distinguish perched water from the regional water table. Use well screens,
  aquifer identity and elevation datums when interpreting observations.
- For constant-density Darcy flow, preserve the definitions
  `h = z + psi`, `q = -K grad(h)` and `u = rho_w g psi`, with documented
  elevation and pressure sign conventions. Saturated Darcy flux is not pore
  velocity. Unsaturated conductivity depends on the hydraulic state.
- For variably saturated flow, document the water-retention and conductivity
  laws, parameter provenance, source/sink terms and storage formulation. Test
  capillary effects, wetting/drying fronts, moving water tables and layered
  contrasts when they control the application. Justify omitted hysteresis,
  preferential flow or lateral unsaturated flow where relevant.
- Standard MODFLOW 6 UZF is a gravity-dominated one-dimensional kinematic-wave
  approximation, not a full Richards-equation solver. DHSVM soil routing also
  must not be described as a general three-dimensional Richards solver. Where
  suction, capillary rise, perched saturation or shallow water-table feedback
  control an outcome, use appropriate physics or demonstrate adequacy against
  a suitable variably saturated reference and observations. An hourly output
  interval does not establish physical adequacy.
- Retain snow mass and energy consistency, precipitation phase and canopy
  effects. Include glacier melt and its routing where material on Rainier.
  Include frozen-soil effects when they materially alter infiltration. Do not
  silently compensate for missing snow/glacier physics by changing aquifer
  properties or recharge multipliers.
- Treat fractured volcanic rock, glacial sediments, confining layers, alluvium
  and hydrothermal systems according to evidence. A first cold-water aquifer
  model does not establish deep volcanic/hydrothermal circulation.
- Record horizontal CRS, vertical datum, length/time units, positive directions
  and flux-versus-volume conventions at every interface. Prefer SI in new
  coupling interfaces; convert legacy feet/day inputs explicitly and test the
  conversions. Use timezone-aware UTC timestamps internally and identify
  accumulation windows; local display times must handle daylight saving time.

## Coupling contract

Before implementing a coupling, document process ownership and exchanges in
`docs/`: model/version, physical domain, grid, time stepping, complete state,
flux units/signs, mapping, feedback, convergence and restart behavior.

- Assign each physical storage volume and flux one owner. DHSVM already
  includes soil moisture and shallow saturated flow. Either retain those in a
  physically distinct shallow domain coupled to a deeper aquifer, or replace
  the overlapping saturated routing. Do not add MODFLOW groundwater underneath
  an overlapping DHSVM aquifer without defining the partition.
- Root-zone drainage is not automatically instantaneous water-table recharge.
  Represent travel and storage through any intervening vadose zone. Place
  the top of a separate vadose module at the actual interface, without
  duplicating soil storage. Return groundwater feedback where it changes
  infiltration, soil saturation, exfiltration or uptake.
- Assign stream routing, gaining/losing exchange, baseflow and ET consistently.
  Count each flux once. Groundwater or vadose ET must consume only the demand
  allocated to that component; do not apply full potential ET repeatedly.
- Map exchanges conservatively using intersecting areas/reaches and integrated
  volumes over the exchange interval. Different grids and time steps are
  allowed. Verify both directions and partially active cells. Interpolation
  of heads alone is not a conservative flux transfer.
- Distinguish numerical time step, coupling interval, output interval,
  observation interval, assimilation cycle and forecast horizon. Substep or
  iterate where rapid responses or strong feedback require it. Demonstrate
  time-step and coupling-interval convergence before adopting asynchronous or
  lagged updates. Hourly interpolation of a daily solution is not an hourly
  groundwater solve.
- In iterative coupling, retain the correct beginning-of-step state; do not
  repeatedly advance storage or reapply a flux while seeking convergence.
  Check head and flux convergence as well as the combined budget.
- Save and restore the complete coupled restart state, including snow/soil,
  vadose memory or wave states, channel/surface storage and aquifer state.
  Demonstrate restart equivalence; restoring heads alone is insufficient.
- When prescribed material parameters change on a fixed grid, retain each
  cell's old conserved `W` and invert the new constitutive law. Do not hold
  pressure fixed and silently change water, clip invalid storage, or modify
  geometry through this transaction. Reject an update outside the declared
  pressure/physics envelope. Recovery restarts must retain separate relaxation
  modes/timescales, event IDs/order, forcing position and the global water ledger.
- Model DataHub groups require exact grid/reference/active-mask identity,
  explicit units, UTC timestamps, valid material epochs/fingerprints and
  finite active physical fields with NaN outside the material domain. Synthetic
  local grids must not claim an earth CRS. Do not run `gaia align` on model
  prognostic/material/event/checkpoint groups: its current horizontal
  reprojection and temporal averaging do not implement conservative 3D coupling.
- Groundwater boundaries need hydrogeologic justification and may cross
  surface watershed divides. Rainier park boundaries, the volcanic edifice,
  a rectangular GIS extent and gauged watersheds are different domains.

## Hazard relevance and resolution

For every application, identify whether each output is an antecedent state,
susceptibility indicator, triggering probability/factor of safety, deformation,
or runout/inundation estimate. Validate the claimed quantity. The time scales
below guide model selection; establish actual steps and meshes through local
process scales, observations and convergence tests.

| Application | Hydrologic role and additional required physics | Scale requirement |
| --- | --- | --- |
| Earthquake liquefaction | Supply antecedent saturation, shallow water-table and pore-pressure profiles. Add shaking demand and geotechnical susceptibility/resistance, using an identified empirical triggering method or an effective-stress model with cyclic excess pore-pressure generation and dissipation. | Hourly antecedent hydrology may initialize a seconds-to-minutes earthquake analysis; resolve susceptible soil layers and site conditions. |
| Lateral spreading | Supply hydrologic conditions to a liquefaction/deformation assessment with stratigraphy, strength, slope and river-bank/free-face geometry. A liquefaction indicator alone does not predict displacement. | Event-scale mechanics and local geometry are required for displacement; regional grids support screening. |
| Rainfall/snowmelt landslides | Supply depth-resolved suction and pore pressure at potential slip surfaces to a stability model with soil/rock strength, slope geometry and relevant root reinforcement. Shallow infinite-slope methods do not establish deep-landslide behavior. Earthquake-triggered failures also need shaking. | Antecedent wetting spans hours to seasons; rapid infiltration may require subhourly steps and hillslope/site refinement. Triggering is distinct from runout. |
| Floods | Supply runoff and streamflow. Depth, velocity, arrival time and inundation require appropriate channel/floodplain hydraulics, including material backwater, structures, storage and boundary conditions. | Test storm and routing steps against catchment response times. Flash flooding can require minutes and local terrain/structure resolution. |
| Wildfire conditioning | Soil moisture and groundwater can inform vegetation water availability. Live/dead fuel moisture, ignition and fire spread require suitable fuel, weather and fire models. | Match the relevant fuel/weather response and spatial support; soil water content is not a direct fuel-moisture measurement. |
| Postfire floods/debris flows | Represent evidence-based changes and recovery in canopy, litter, infiltration, water repellency, roughness, erodibility and roots. Distinguish runoff-driven erosion/debris flows from infiltration-driven slope failure; use appropriate initiation and runout models. | Short rainfall bursts, often evaluated over 5–30 minutes, can control initiation. Hourly means alone cannot establish peak-intensity hazard. |

Do not infer site-scale hazard skill merely by downscaling a coarse hydrologic
field onto a fine DEM. Regional screening, watershed forecasting and
site-specific engineering analysis require different evidence. Do not assume
burning always decreases hydraulic conductivity; require local evidence or
test competing plausible scenarios. Account for compound events where material
(rain-on-snow, wet antecedent conditions plus shaking, and postfire storms).

## Observations and data assimilation

- Inventory observation location, sensor depth or well screen, units, datum,
  uncertainty, spatial support, quality flags, sampling frequency, arrival
  latency and actual uptime. Preserve observation time and arrival time.
- Use observation operators appropriate to sensor support: near-surface
  moisture is not a full profile; a confined well is not the shallow water
  table; discharge integrates upstream processes. Do not infer deep aquifer
  structure uniquely from surface moisture and outlet discharge.
- Treat precipitation primarily as uncertain forcing. Retain gauge undercatch,
  rain/snow phase, elevation and radar uncertainty where relevant. Snow-energy
  modeling also requires temperature, humidity, wind and radiation or justified
  estimates, not precipitation alone. Seek snow/SWE and depth-profile or
  piezometer observations where needed to resolve ambiguity.
- Avoid treating observations already used in forcing products as independent
  evidence without accounting for shared errors. Do not silently replace
  missing precipitation, moisture or head with zero.
- Separate offline parameter estimation and spin-up from online state updates.
  Use ensembles, localization and physically justified parameter uncertainty;
  demonstrate ensemble-size and localization sensitivity. Avoid dense
  state-by-state covariance matrices for regional applications.
- Enforce constitutive consistency and physical bounds after assimilation.
  Record all analysis-induced storage changes as explicit increments in the
  water ledger. Clipping an invalid state is not a water-conservation method.
  Reconcile pressure, moisture and storage consistently before advancing.
- Support delayed observations through a documented replay or fixed-lag
  approach. Hourly sampling does not guarantee hourly telemetry. Products must
  expose their valid time, issue time, observation freshness, missing-data
  status and uncertainty.
- Hold out sensors, sites and events. Report forecast skill as well as fit at
  assimilated observations, event peaks/timing, head/moisture errors and
  uncertainty calibration. Hazard probabilities also need reliability,
  false-alarm and missed-event evaluation against independent events.

## Performance and core-code optimization

Core DHSVM and MODFLOW changes are within project scope when they improve
coupling, runtime or memory while preserving the required physics. Prefer
existing supported interfaces when adequate; source modifications are allowed
when justified. Keep upstream versions/commits and project patches traceable.

- Profile representative cases before optimizing. Measure wall time, peak
  resident memory, I/O, solver/coupling iterations and ensemble throughput on
  the actual workstation. Separate warm-up/build time from recurring cycles.
- Begin with persistent model instances, in-memory exchange through documented
  APIs, precomputed conservative mappings, buffer reuse, selective outputs,
  efficient forcing access and ensemble batching. DHSVM may require a
  reentrant step/state interface; MODFLOW 6 already provides an API. Do not
  restart whole historical simulations or write full spatial dumps every hour
  merely to exchange small arrays.
- Verify reentrancy and ownership of mutable state before parallel execution.
  Prefer isolated ensemble processes when global/static state is unsafe.
  Avoid nested parallel oversubscription. Test initialization, repeated runs,
  checkpoint/restart and cleanup, including failure recovery.
- Verify API behavior against the pinned release, including allowed state
  writes, pointer lifetimes, reinitialization and time-step control. A function
  argument does not guarantee that runtime time-step changes are honored.
  Synchronize actual model time with forcing and assimilation windows.
- Cache only demonstrably invariant quantities. Do not freeze head-dependent
  conductances, saturation-dependent storage or dynamic routing geometry.
  Parallel routing must avoid races when multiple cells update a shared
  neighbor or channel; verify local and global conservation after refactoring.
- Evaluate allocation/layout improvements, vectorization, sparse solver and
  preconditioner choices, domain decomposition and selective refinement based
  on profiles. GPU acceleration, reduced precision, surrogate models and
  longer deep-groundwater steps require demonstrated benefit and validation;
  they are not assumed speedups.
- Compare every meaningful optimization with a pinned reference on identical
  inputs/hardware. Check water/energy budgets, solver convergence, field/state
  differences and event-relevant outputs, including pressure near slip
  surfaces, water-table crossings and flood peaks. Accept defined numerical
  equivalence rather than requiring bitwise equality across architectures.
- Do not achieve speed by silently weakening physical processes, convergence
  tolerances, calibration, ensemble uncertainty or hazard resolution. Document
  any approximation, its error and its validated operating range.
- Benchmark storm onset, dry-to-wet transitions, rain-on-snow, near-surface
  groundwater and representative postfire cases, not only benign averages.
  Report tail cycle times, solver failures and recovery behavior.
- Distinguish one-hour nowcast propagation from multi-day forecast cost.
  Include ingestion, observation processing, coupling, assimilation, output,
  restart and any delayed-data replay in operational timing. Projections are
  not measured performance; HPC results do not establish workstation capacity.

## Development and acceptance checks

- Keep raw/reference data immutable and record download dates, versions,
  licenses, transformations, configuration and model/software commits.
- Preserve the distinction between the published calibrated MODFLOW-NWT
  reference and the existing MODFLOW 6 template. Changing engines or temporal
  resolution requires renewed physical and numerical validation.
- The current lightweight test command is
  `PYTHONPATH=src python3 -m unittest discover -s tests`; `pixi run test` uses
  the project environment. MODFLOW smoke runs require the relevant executable.
- Gaia interoperability is an optional `datahub` dependency extra pinned in
  `pyproject.toml`; it is intentionally outside the default numerical Pixi
  environment. The extra requires Python 3.12+ because the pinned Gaia stack's
  `rioxarray>=0.23` requires it; the independent core supports Python 3.11+.
  Its CI job sets `HYDROSTATE_REQUIRE_GAIA=1` and must exercise
  actual Gaia write/read APIs, including complete recovery-memory checkpoints.
  Generic DataHub timestamps preserve nanoseconds; the present hydraulic and
  recovery integration clocks support microseconds and must reject finer
  precision rather than silently truncating it.
- For scientific/coupling changes, add meaningful benchmarks as applicable:
  hydrostatic/no-flow equilibrium; infiltration/drainage against a suitable
  reference; layered/perched and moving-water-table behavior; gaining/losing
  streams; conservative mapping; restart equivalence; and assimilation budget
  closure. Test mesh, vertical layering, time-step and coupling sensitivity.
- Software unit tests, physical benchmarks, calibration and independent
  validation are different checks. Report which were performed. Documentation
  changes do not require new numerical tests.
- Before a scientific or hazard claim, document domain boundaries, physical
  adequacy, parameter provenance, calibration/validation, resolution sensitivity,
  uncertainty, budget closure and measured operational performance appropriate
  to that claim. Present limitations directly; do not imply operational warning
  capability from an unvalidated prototype.

## Primary references

Use primary literature and official documentation; verify current APIs and
research-fork status before implementation. These are starting points:

- [PNNL DHSVM source](https://github.com/pnnl/DHSVM-PNNL) and
  [documentation](https://www.pnnl.gov/projects/distributed-hydrology-soil-vegetation-model/documentation).
- [USGS MODFLOW 6 groundwater-flow theory](https://pubs.usgs.gov/tm/06/a55/tm6a55.pdf),
  [UZF input reference](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-uzf.html)
  and [MODFLOW API](https://doi.org/10.1016/j.envsoft.2021.105257).
- [DHSVM dynamic stream/groundwater development](https://doi.org/10.5194/hess-30-5145-2026).
- [Mass-constrained soil-moisture assimilation](https://doi.org/10.5194/hess-16-105-2012).
- [USGS liquefaction mechanisms and lateral spreading](https://www.usgs.gov/programs/earthquake-hazards/what-are-effects-earthquakes).
- [USGS TRIGRS infiltration/slope stability](https://pubs.usgs.gov/of/2008/1159/).
- [USGS postfire short-duration rainfall analysis](https://pubs.usgs.gov/publication/ofr20231025/full).
- [USGS wildfire effects on Oregon Cascades soil hydraulics](https://www.usgs.gov/publications/effects-wildfire-soil-hydraulic-properties-western-oregon-cascades).
- Regional hydrogeology and reference-model sources are also recorded in
  `docs/data_sources.md`; this is a starting inventory, not proof of suitability
  at hourly or hazard scales.
