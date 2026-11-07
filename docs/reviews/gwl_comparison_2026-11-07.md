# Review of the regional soil twin and the Hydrostate–Mechanicalstate approach

Review date: 7 November 2026. This is a code-and-equations comparison and a
revision of the proposed architecture. It does not certify a regional forecast,
a hazard adapter or the reported field calibration. No upstream GWL source or
observational/cloud data was modified.

## Reviewed sources and scope

| Repository | Reviewed physical-code revision |
| --- | --- |
| [GAIA Digital Twin of Soil / gwl-space-time-smooth](https://github.com/gaia-hazlab/gwl-space-time-smooth) | `435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f` |
| [Hydrostate](https://github.com/gaia-hazlab/hydrostate) | `c9328f331c6e6985054c0651bde15cc72235fb64` |
| [Mechanicalstate](https://github.com/gaia-hazlab/mechanicalstate) | `b86d3b2e8a138f238e9db94dd7ece1c6c89db3eb` |

The groundwater repository is **not merely a spatial smoother**. Its current
source includes a reduced daily soil–groundwater budget, snow and petrophysical
functions, hysteresis, footprint-aware snapshot analysis, uncertainty/resolution
diagnostics and forecast/export workflows. The roadmap reports a 2.96-million-cell
western Cascades grid at 90 m, field-derived calibration and explicitly multiscale
hazard consumers. Those published numerical results were inspected, not independently
regenerated here. Older monthly/1 km README and instruction descriptions are not a
reliable substitute for actual current source and versioned results.

We inspected source and tests, ran 31 selected existing offline tests, and
independently reproduced a combined-budget counterexample. We did not rerun the
full geospatial suite, field calibration, large-data pipeline, FuXi GPU adapter or
actual downstream hazard consumers. Official Landlab source/API and primary hazard
references were also inspected. A format exporter is distinguished from execution
of a scientific consumer.

## Comparison

| Topic | Regional soil twin | Hydrostate + Mechanicalstate | Revised implication |
| --- | --- | --- | --- |
| Main strength | Observation/data workflows, regional priors, native support and analysis diagnostics | Explicit local variably saturated flow, conservative material transactions and complete damage/recovery checkpoints | Reuse regional orchestration; offer selectable process engines |
| Hydrology | Reduced root bucket + groundwater reservoir, slope partition, HAND/TOPMODEL diagnostics | Cartesian 3D equilibrium-retention Richards prototype with conservative face fluxes | Compare within regime; do not equate a reduced state with a resolved pressure profile |
| Groundwater analysis | Snapshot BLUE on phreatic anomaly, with screened-head semantics and residual/anomaly Matérn prior | No ensemble/filter analysis currently implemented | Reuse and extend the existing analysis machinery rather than duplicate it |
| Observation support | Existing footprint operators, seismic sensitivity exploration, resolution/DFS machinery | Separate direct-sensor CSV comparisons; acquisition-specific dv/v/GNSS-IR operators remain incomplete | Preserve native supports and shared errors in the common inference layer |
| Memory and events | Hydraulic/seismic hysteresis work exists; full cycling/damage integration remains incomplete | Independent prescribed hydraulic/elastic anomalies, finite recovery modes and exact UTC events | Assign hydraulic constitutive memory one owner; share material/recovery capability |
| Forecast restart | Reduced forecasts reset important initial states in current calls | Actual hydraulic/material/schedule restart and Gaia continuation tested on a tiny synthetic domain | Extend full state restart to regional engines, including snow/analysis state |
| Probability | Posterior marginal variances in the matrix-free analysis path; no joint posterior samples there | No calibrated joint posterior; one synthetic material scenario | Neither currently supplies validated coherent joint hazard ensembles |
| Delivery | 2D regional grids, source/support metadata, Landlab ASCII/COG exporter | Strict 3D same-grid Gaia state/material/checkpoint extension | Add a reduced-product/boundary protocol; keep restart identity strict |
| Hazards | Explicit target consumers and export plans; mechanics scaffold | Antecedent state/material diagnostics; no strength/deformation/hazard solver | Name and execute each actual consumer before claiming coupling |
| Operation | Larger real-data workflow and performance work; forecast/GPU and dense observation-side costs remain | Tiny synthetic runtime and verification; regional throughput unmeasured | Measure end-to-end ensemble cycle and tail latency on the workstation |

## Findings that change the design

### 1. Confirmed combined capillary-transfer imbalance

In the reviewed [water budget](https://github.com/gaia-hazlab/gwl-space-time-smooth/blob/435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f/src/models/water_budget.py#L214),
capillary rise is credited to root storage. The groundwater update at lines
266–278 includes recharge and drainage but no matching capillary debit.
For root storage `S` in mm and aquifer head anomaly `h` in m, the intended
combined balance is

```text
Delta(S + 1000*Sy*h) = P - ET - runoff - interflow - baseflow.
```

The current implemented right-hand side acquires an additional positive capillary
term. A one-cell case with zero P/PET, 1 m root depth, theta_wp=.10,
theta_fc=.28, theta_sat=.45, Sy=.30 and initial DTW=.50 m, initialized at wilting
point for one day, produces:

| Quantity | Executed result |
| --- | ---: |
| Soil storage change | +0.739233196 mm |
| Aquifer storage change | 0 mm |
| External inputs and outputs | 0 mm |
| Recorded internal capillary transfer | +0.739234865 mm |
| Combined residual, expected zero | +0.739233196 mm |

The small rounding difference comes from float32 output. Analytic prescribed
capillary supply is `30/30.436875 * .75 = .739234892 mm`. A closed internal transfer
must debit the aquifer as well as credit the soil. This is a confirmed defect
for the reproduced case; it does not alone establish the magnitude of errors in
other cases or invalidate every reported field comparison.

Reproduce from a trusted checkout at the pinned revision:

```bash
python scripts/review_gwl_capillary.py \
  --checkout /path/to/gwl-space-time-smooth \
  --output runs/reviews/gwl_capillary_counterexample.json
```

The [committed diagnostic](gwl_capillary_counterexample.json) records source hash,
inputs, tolerance and `combined_budget_passed=false`. Exit zero from the script
means the audit executed, not that the model passed conservation. Existing
`test_water_budget.py` root-zone checks do not test this combined closed transfer.
A second [non-overlapping-column counterexample](gwl_capillary_nonoverlap_counterexample.json)
uses DTW=1.50 m, below the 1 m root-zone base. It still creates `0.246414542 mm`
with no external water, so the finding is not dependent on root/aquifer geometric
overlap in the shallower case. Run with `--initial-dtw-m 1.5` to reproduce it.
A repair must check storage representability, shallow-table feedback and renewed
calibration, rather than merely subtracting an arbitrary clipped amount.

### 2. Forecasts do not yet continue the complete posterior

The reviewed [forecast](https://github.com/gaia-hazlab/gwl-space-time-smooth/blob/435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f/src/models/forecast.py#L105)
initializes snow at zero (line 123); the budget initializes root storage from
field capacity/wilting point and head anomaly at zero (water_budget.py155–158).
Hysteresis is initialized from the beginning of a supplied trajectory. These
functions lack complete input/output restart of accepted posterior stores and
memory. `theta_ref=wb.theta[0]` at forecast.py162 is after the first hydrological
step, although comments describe an initial analysis reference.

The repository itself correctly distinguishes current snapshot BLUE from target
storage-first recursion in
[hydro_state.py49–61](https://github.com/gaia-hazlab/gwl-space-time-smooth/blob/435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f/src/models/hydro_state.py#L49).
A forecast cannot be called posterior cycling until it resumes the complete
accepted state, including pre-step snow and seismic/constitutive reference.
Do not substitute the posterior head for the fixed HAND drainage-reference
parameter without separating those physical roles.

### 3. Preserve joint uncertainty and the actual likelihood

[observability.py151–153](https://github.com/gaia-hazlab/gwl-space-time-smooth/blob/435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f/src/models/observability.py#L151)
explicitly says the matrix-free path supplies posterior marginal variances, not
posterior samples or cross-cell posterior covariance. Independent draws from
mean/sigma cannot supply the joint distribution needed by routed or spatially
compound hazards. Stable member identities must span water, materials, damage,
forcing, lead time and hazard evaluation.

BLUE constructs diagonal observation error from a variance vector; the seismic
band inversion also uses a diagonal error representation and its own prior.
Shared seismic references and retrieval priors can therefore be recycled unless
explicitly handled. Prefer forward prediction of measured bands from the joint
state, or retain retrieval prior, averaging kernel and cross-covariance. An
atmospheric GNSS product and GNSS-IR shallow moisture need separate operators.

The exact FFT prior and support/resolution diagnostics are important reusable
assets. They do not eliminate the dense `n_cell x n_obs` intermediate documented
in observability.py140–144. Full-domain analysis/ensemble throughput needs direct
measurement; a fast prior multiplication is not a complete cycle benchmark.

### 4. Reduced regional hydrology and 3D physics have different roles

The regional budget partitions root drainage into instantaneous recharge and
interflow. It does not resolve neighboring-cell Darcy exchange or intermediate
vadose transit/storage in that routine. TOPMODEL is an equilibrium diagnostic;
`accumulate_runoff` is explicitly a weighted proxy rather than conservative
network routing. These approximations can support evaluated regional antecedent
products, but they do not identify transient pressure on a perched slip surface.

Hydrostate supplies a useful reference for such local regimes, but our current
solver has its own large gaps: synthetic grid/time evidence, no established
nonlinear solution mesh convergence, isotropic toy materials, no complete
DHSVM/MODFLOW integration, no operational ensemble analysis, and no regional
performance result. Mechanicalstate's constant illustrative net stress and
granular frame are not general stress/strength solutions. Damage and fire
amplitudes/timescales remain prescribed; geometry-changing events are rejected.
A new universal statewide 3D solver is not justified by either implementation.

### 5. Equal containers are not equal physical states

Regional products include DTW/WTE and uncertainty on 2D grids. The existing
Hydrostate–Mechanicalstate protocol carries full 3D hydraulic/material arrays,
clock/epoch identity and complete owner checkpoint. A regional raster is a
prior/observation/boundary candidate, not a 3D restart.

Specific hazards at the interface include `mm per step` versus `mm/day`,
pressure-based VG alpha (`1/kPa`) versus head-based alpha (`1/m`), mixed datum or
screened-head semantics, and coarse temporal forcing statistically delivered at
90 m. The exporter uses bilinear reprojection for fields including flux rates
and sigmas; neither water conservation nor covariance propagation follows from
that operation. Its missing-CRS fallback should not be treated as georeferencing
evidence. Derive actual grid spacing and preserve native support; metadata
currently also contains legacy resolution defaults.

### 6. Named hazard consumers need equation-level integration tests

Current official [Landlab API/source](https://landlab.csdms.io/generated/api/landlab.components.landslides.landslide_probability.html)
shows `LandslideProbability` deriving relative wetness from recharge,
transmissivity and topography. Its ordinary fields do not consume externally
provided transient DTW or slip-surface pore pressure. The original study uses
annual maximum daily recharge with steady-state subsurface flow to estimate
annual initiation probability:
[Strauch et al. (2018)](https://doi.org/10.5194/esurf-6-49-2018).
The GWL exporter tests parse ASCII rather than instantiate the consumer. Keep
that regional screening interface separate from transient pressure-driven
stability, and demonstrate the actual chosen component consuming the fields.

The [Sanger/Geyin/Maurer preprint](https://arxiv.org/abs/2509.10962), especially
pp. 29–30 of the reviewed PDF, describes published premapped coefficients
used with magnitude-scaled PGA and identifies groundwater depth as a future
user-input extension. Baseline groundwater affecting map construction does
not make current DTW an input to static map evaluation. A dynamic update needs
a demonstrated compatible surrogate/rerun or another identified geotechnical
model. Susceptibility, conditional triggering and lateral-spread displacement
remain separate quantities.

Flood source volumes must enter a reach/network routing and hydraulic model;
weighted runoff or groundwater boundary outflow is not an observed hydrograph.
The [updated architecture](../soil_state_architecture.md) lists each hazard's
required pressure/depth, forcing, strength, geometry, temporal and spatial support.

### 7. Independent validation must respect issue time and calibration dependence

`gwl_dynamic.py94–97` centers wells using the whole requested time window. This
is suitable for retrospective analysis; an operational hindcast must freeze
that baseline, data availability and all calibration at issue time. Similar
care is needed for seismic references and observation latency.

Gauge BFI, Q/P and seasonal head amplitude provide valuable distinct constraints,
but fitting different parameters to different metrics does not prove statistical
independence, uniqueness or held-out forecast skill. Require sensitivity,
equifinality and spatial/event holdouts. Surface-versus-deep moisture, confined
versus phreatic head and routed versus generated flux must be scored separately.

## Executed checks and proposed promotion tests

Existing selected suite run in the checkout:

```bash
PYTHONPATH=. /path/to/hydrostate-python -m pytest -q \
  tests/test_water_budget.py tests/test_forecast.py \
  tests/test_hydro_state.py tests/test_hysteresis.py
```

The review recorded 31 passing tests. The additional independent closed-transfer
counterexample fails combined conservation. A passing existing suite therefore
does not settle that missing invariant. The standalone diagnostic was rerun by
the primary reviewer with the committed source hash and values above.

Priority promotion tests are complete restart/replay through snow, wetting/drying,
disturbance and analysis; constrained analysis-water increments; mixed-grid/time
volume conservation; raw-versus-derived observation error accounting; coherent
hazard ensembles; actual consumer calls; withheld-issue/event forecasts; and
measured complete workstation cycles. See the
[revised decision and gates](../soil_state_architecture.md) and independent review
prompts for [Fable](fable_adversarial_prompt.md) and
[Antigravity](antigravity_adversarial_prompt.md).
