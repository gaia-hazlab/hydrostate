# Shared soil state for regional analysis and hazard forecasts

Decision date: 7 November 2026. This updates the architecture after reviewing
`gwl-space-time-smooth`, now titled **GAIA Digital Twin of Soil**. It supersedes
the suggestion to build an independent `SoilStateRunner` beside that project.
It is an architecture decision and acceptance plan, not a newly implemented
assimilation service, nested model or hazard forecast.

The [pinned comparison](reviews/gwl_comparison_2026-11-07.md) records source,
executed checks and unresolved findings. Independent review prompts are available
for [Fable](reviews/fable_adversarial_prompt.md) and
[Antigravity](reviews/antigravity_adversarial_prompt.md).

## Decision and division of responsibility

Use a **shared, model-agnostic soil-state protocol**, with the existing GAIA soil
twin as the proposed regional analysis/orchestration host. Hydrostate and
Mechanicalstate provide selectable physical process capabilities and verification.
A shared 90 m reporting grid can be useful; it is not a universal process mesh,
observation resolution or hazard resolution. These roles are an integration
proposal; the three projects are not yet connected end to end.

| Owner | Responsibility | Present boundary |
| --- | --- | --- |
| GAIA soil twin / `gwl-space-time-smooth` | Regional observation catalog, footprint operators, geological/statistical priors, analysis diagnostics, forcing/data workflows and intended nowcast/forecast/hazard orchestration | Existing snapshot analyses and reduced forecasts do not yet constitute complete posterior cycling; confirmed budget/restart debts must be addressed |
| Hydrostate | Water stores, pressure, constitutive consistency and conservative fluxes in its assigned domain; reference Richards tests and optional refined watershed/site engine | Synthetic independent Richards prototype; full DHSVM–MODFLOW driver, regional performance and nonlinear mesh adequacy remain unverified |
| Mechanicalstate | Structural priors, physical porosity, material laws, independent disturbance/recovery memory and elastic observation predictions | Prescribed fixed-grid scenarios; no calibrated damage-generation, strength or deformation solver |
| Gaia data plane | Version observations, supports, units, state products, uncertainty structures, source lineage and complete owner checkpoints | Shared DataTree/Zarr serialization is not a physical, statistical or hazard interface by itself |
| Named hazard consumers | Stability/triggering, earthquake response, routing/hydraulics, displacement, runout or fire processes | Each consumer needs its own actual API, physical assumptions, scales and validation |

Do not duplicate the soil twin's observation/inference machinery in a separate
regional controller. Reuse it after testing its operators and statistical
assumptions; place hydraulic backends behind a documented advance/analyze/
checkpoint interface. Avoid making a new regional 3D solver the prerequisite
for every useful observation-anchored soil product.

## Common state, several physical representations

A complete accepted member contains all prognostic stores, constitutive memory,
material/recovery modes, forcing position, references and accumulated ledgers
owned by its active engines. Derived moisture, pressure, saturation and elastic
fields must agree with that state. Do not independently assimilate inconsistent
versions of the same physical degrees of freedom.

The regional reduced model may carry root-zone storage and unconfined aquifer
storage. A Richards patch carries pressure and the conserved coordinate
`W=theta(psi)+Ss*max(psi,0)`. Specific yield in a reduced reservoir is not an
extra storage term to add to the same water-table drainage already represented
by Richards retention. Snow, surface, channel and deep-aquifer stores each have
one owner, including when DHSVM or MODFLOW is introduced.

Distinguish:

- Phreatic water-table elevation, positive-down water-table depth, screened
  aquifer head, perched saturation and pressure at a specified depth.
- Physical liquid content `theta`, physical saturation `theta/phi`, retention
  effective saturation and the elastic contribution to conserved `W`.
- Equilibrium retention, hydraulic scanning-curve hysteresis, rate-dependent
  capillarity and persistent structural damage. A diagnostic seismic memory
  law does not automatically change hydraulic retention or flux.
- Small-strain stiffness, geotechnical resistance, triggering and deformation.
  A modulus-loss variable does not supply CPT resistance or residual strength.

A regional DTW raster does not uniquely determine a 3D pressure/moisture profile.
Initializing a refined patch therefore needs an explicit lifting/conditioning
procedure using stratigraphy, boundary heads/fluxes, profile observations,
spin-up and uncertainty. Name hydrostatic assumptions when used. In perched,
confined, strongly transient or sloping flow, reject unsupported pressure
reconstruction instead of silently applying a hydrostatic profile.

## Two interfaces, kept separate

The existing Hydrostate–Mechanicalstate fixed-grid checkpoint protocol remains
strict. Exact coordinates/reference/mask/material identity protects a running
accepted model. It should not be loosened to accept statistically downscaled
regional rasters as solver restarts.

Add a separate **reduced-product / boundary / observation envelope** for regional
products and hazard handoffs. A future executable validator should require:

| Information | Required meaning |
| --- | --- |
| Quantity and role | Named DTW, phreatic elevation, screened head, pressure, moisture, flux or volume; observation/prior/boundary/forecast role |
| Units and signs | SI internally; explicit positive direction; rates versus per-interval totals |
| Reference and domain | Horizontal CRS, vertical datum, aquifer/screen identity, surface definition and physical active/coverage masks |
| Clocks | Valid time, issue time, observation arrival, forecast lead, interval bounds and averaging/integration convention |
| Support | Native spatial, depth and temporal sensitivity; delivered grid; transformations, coverage and missing values |
| Probability object | Deterministic field, marginal summary, covariance-bearing distribution or coherent ensemble; no implicit conversion between them |
| Ensemble identity | Stable ensemble/member identifiers aligned across water, materials, forcing, damage, times and hazard calculations |
| Lineage and quality | Source/code revisions, data used in calibration/assimilation, confidence, applicability and freshness |
| Consumer declaration | Actual model/version/API, required fields/support, assumptions and supported output category |

Supported conversions must be explicit: DTW to phreatic elevation needs same-datum
surface elevation; mm per interval to m/s needs exact interval duration; pressure
head to pressure needs declared fluid density and gravity. GWL's pressure-based
VG alpha in `1/kPa` converts to Hydrostate's head-based `1/m` through
`alpha_head=alpha_pressure*rho_w*g/1000`; direct parameter copying is invalid.
GNSS atmospheric ZTD/TEC, GNSS-IR moisture and GNSS displacement require distinct
measurement operators.

Do not use bilinear interpolation as a conservative water-volume transfer.
For flux coupling, integrate overlap areas and actual time bounds, handle partial
coverage explicitly, and verify source/receiver volumes. Pressure/head mapping
must respect the governing problem and hydrostratigraphy. Fine-grid material
conditioning requires a validated model and uncertainty; it does not create
observed fine-scale hydrological information.

## Analysis, forecast and probability

Use one complete posterior and one issue-time definition. Regional and local
engines can have different meshes and steps, but must share accepted exchange
boundaries, event timing, units and ledger ownership.

1. Ingest observations at native supports and actual windows, with arrival time.
2. Forecast each coherent member using uncertain forcing and material states.
3. Predict raw measurements through appropriate forward operators; preserve
   seismic reference, frequency, path/depth kernels and cross-errors.
4. Analyze with an evaluated filter/smoother; select localization and observation
   batching from sensitivity/observability and measured runtime.
5. Reconcile constitutive state and bounds. Material-only changes preserve old
   water. Observation analysis records
   `DeltaV_analysis=sum(V_i*(W_i_analysis-W_i_forecast))`; changes to other
   stores have corresponding explicit entries. Clipping is not a balance.
6. Checkpoint every owned store/memory, future event, random state/member ID,
   forcing position, parameter field and cumulative budget.
7. Propagate to hazard-specific consumers and score withheld observations/events.

Do not assimilate an already analyzed groundwater or seismic retrieval as an
independent observation of the raw data that produced it. Either use original
measurements, preserve retrieval prior/averaging kernels/covariance, or adopt a
single clearly named posterior. Freeze priors, references, preprocessing and
calibration at forecast issue time to prevent future leakage.

A marginal mean and sigma field does not determine a joint distribution.
Independent pixel/time draws invent spatial, temporal and water–damage–strength
dependence. Hazard probabilities require coherent samples or another explicitly
specified, justified uncertainty representation. Preserve identical member
identity through all downstream calculations and shared event drivers.

## Hazard-specific products and scales

The following are selection and verification requirements, not universal grid or
step thresholds. Output interval, internal step, sensor cadence, forecast issue
cadence and useful lead time are different quantities.

| Application | Required handoff and added physics | Scale and acceptance |
| --- | --- | --- |
| Regional shallow-landslide screening | Recharge distribution plus compatible transmissivity, soil thickness, slope, cohesion/root strength and friction for the named Landlab component | Its recharge-derived wetness and steady-state assumptions must be retained and validated; an annual/climatological initiation probability is not automatically event-time probability |
| Storm-time shallow landslides | Transient pressure/suction at candidate slip depths, density and root/material strength; appropriate stability component | Hours–season antecedent states; shorter steps where wetting fronts matter; hillslope and vertical-layer refinement must resolve geometry and perched pressure |
| Deep landslide/runout | Relevant deep/fracture pressure and geology, actual slip geometry, strength/deformation and runout solver | Site evidence and identified processes; shallow infinite-slope capability is insufficient |
| Regional liquefaction screening | Appropriate shallow unconfined DTW distribution, geology/material proxies and a compatible shaking-conditioned model | Regional susceptibility differs from event triggering; dynamically changing DTW must enter an actual validated model input or documented coefficient rerun |
| Site liquefaction / lateral spreading | Saturation and pore-pressure/stress profiles in susceptible layers, CPT/SPT or justified Vs resistance, shaking demand, bank/free-face geometry and displacement model | Hourly antecedent hydrology may initialize seconds-to-minutes earthquake analysis; resolve susceptible layers and validate the claimed triggering/displacement quantity |
| River floods/inundation | Conservatively mapped runoff/interflow/baseflow source volumes/timing and owned channel states; reach routing and channel/floodplain hydraulics | Network/reach and floodplain meshes; internal hydraulic steps selected through stability/accuracy and space/time convergence, not the soil reporting grid |
| Wildfire/postfire | Root-zone water availability, independently constrained fuel/surface/root properties and recovery; rainfall burst, erosion/debris-flow or fire model as needed | Hourly means can erase controlling 15-minute rainfall intensity; soil moisture is not fuel moisture and hydraulic recovery is not ecological recovery |

Two downstream assumptions specifically require correction before promotion:

- Current Landlab `LandslideProbability` reconstructs wetness from recharge,
  transmissivity and topography; its ordinary API does not consume external DTW
  or transient slip-plane pressure. Keep recharge-driven screening separate
  from a pressure-driven transient stability adapter. Run the actual consumer
  in integration tests, rather than only parsing its exported ASCII files.
- Published Sanger/Geyin/Maurer premapped liquefaction coefficients are evaluated
  with magnitude-scaled PGA. Groundwater depth affects the premapped model,
  but a dynamic groundwater input is not automatically available in that map
  interface. Verify the exact implementation/version; use a validated update
  pathway or another appropriate geotechnical model before calling this a
  dynamic groundwater-conditioned forecast.

Primary references:
[Landlab API](https://landlab.csdms.io/generated/api/landlab.components.landslides.landslide_probability.html),
[Strauch et al. (2018)](https://doi.org/10.5194/esurf-6-49-2018),
[Sanger/Geyin/Maurer model preprint](https://arxiv.org/abs/2509.10962),
[Caltrans liquefaction guidance](https://dot.ca.gov/-/media/dot-media/programs/engineering/documents/geotechnical-services/202001-gm-liquefactionevaluation-a11y.pdf),
[HEC-RAS grid/time guidance](https://www.hec.usace.army.mil/confluence/rasdocs/r2dum/latest/running-a-model-with-2d-flow-areas/selecting-an-appropriate-grid-size-and-time-step),
[USGS postfire report](https://pubs.usgs.gov/publication/ofr20231025/full).

## Sequence and gates

1. **Close identified physics/cycling gaps.** Reproduce and fix combined
   soil–aquifer capillary accounting; define complete reduced-model restart,
   snow/hysteresis memory and immutable reference state. Keep legacy calibration
   separate from repaired daily/subdaily physics and revalidate it.
2. **One executable vertical slice.** One catchment and a named hazard consumer,
   direct sensors, bounded uncertain parameters, complete analysis/forecast,
   observed versus forecast issue-time discipline and explicit DA water ledger.
3. **Add seismic/GNSS information.** Validate native operators and shared errors;
   compare with direct-sensor-only forecasts using withheld data, alternative
   truth physics and identifiability tests.
4. **Refined process patches where justified.** Compare reduced and Richards
   engines on layered/perched/water-table cases, then introduce conservative
   boundary exchange, constitutive memory and exact restart/replay tests.
5. **Hazard probability and operational promotion.** Test joint ensemble
   dependence, threshold/arrival/peak errors and reliability with independent
   events. Benchmark ingestion, analysis, ensemble propagation, replay, mapping,
   hazard calls and I/O on the actual workstation, including peak memory and
   tail latency. Regional operation remains a measured acceptance condition.

Acceptance requires local and combined budgets, continuous-versus-restarted
snow/wetting/drying/disturbance/analysis equivalence, conservative mixed-grid
and irregular-time exchanges, actual consumer execution, and field forecast
skill appropriate to each hazard. A pass on schema, a small water residual,
a calibrated seasonal amplitude or a 70-cell runtime proves only its stated
scope. Budget closure and high resolution labels do not certify hazard skill.
