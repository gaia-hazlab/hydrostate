# Reproducible problems and use cases

These experiments progress from known mathematical solutions to synthetic
watershed processes and sensor comparisons. The 3D valley is a new reference
experiment; the older 2D hillslope contains additional surface/ecological
processes. Their different scopes are deliberate and recorded below.

## 1. Verify a known flow before interpreting a watershed

```bash
hydrostate verify --output runs/verification/report.json
python -m unittest discover -s tests -v
```

The [technical reports](technical/) connect governing equations to the actual
finite-volume kernel, analytical derivatives, Darcy fluxes, hydrostatic
equilibrium, conservation and manufactured solutions. Constant-conductivity,
linear-storage tests measure spatial and temporal order on the same face graph,
boundary assembly and time integrator used by the nonlinear solver. Nonlinear
checks have their own stated coverage; passing the linear test does not verify
all nonlinear wetting fronts, water-table crossings or surface boundaries.

## 2. Three-dimensional valley: pressure, moisture and a groundwater boundary

```bash
hydrostate valley --hours 72 --output runs/valley3d/results.json
hydrostate sensors --model runs/valley3d/results.json \
  --sensors examples/sensors/sensors.csv \
  --observations examples/sensors/observations.csv \
  --output runs/sensors/report.json
python scripts/plot_valley.py --output runs/valley3d/valley3d.png
```

All numbers in this example are hypothetical. The domain is 100 × 80 m with
0–6 m elevation and 10 × 8 × 12 Cartesian cells, of which 724 are active.
Horizontal spacing is 10 m; vertical spacing is 0.5 m. The target surface
declines eastward and rises away from the valley centre:

```text
z_target(x,y) = 3.5 + 0.015 (100 - x) + 0.025 |y - 40|  [m]
```

Only complete cells below this surface are retained. This staircase geometry
uses true Cartesian distances and elevation, rather than an uncorrected
terrain-following stencil. Its geometric error still requires refinement.
Soil `Ks = 2e-5 m/s`, lower substrate `3e-6 m/s`, and a finite central lens
`3e-7 m/s`; retention parameters are recorded in the code. These values have
no field calibration. The lens is a conductivity-contrast test, not evidence
of a perched aquifer in this particular run.

Initial total head slopes weakly in both horizontal directions. Fixed total
head `H = 2.3 m` applies only at the outer east boundary, `x = 100 m`.
The bottom, west, north/south and vertical staircase risers are impermeable.
Liquid infiltration is prescribed on horizontal top faces:

| Elapsed interval (h) | Liquid rainfall/infiltration (mm/h) |
| --- | --- |
| 0–12 | 0.05 |
| 12–30 | 0.20 |
| 30–48 | 0.05 |
| 48–60 | 0.15 |
| 60–72 | 0.05 |

At every accepted substep, the code reconstructs the top-face pressure needed
to carry this prescribed flux. It fails if the boundary requires ponding.
This restricted example consequently avoids pretending that all storm rain
can infiltrate. It omits snow, canopy, roots/ET, surface runoff and channel
routing. East-boundary discharge is a groundwater boundary flux; the sensor
example's `streamflow` kind exercises interval/units/outlet matching, not a
validated river-discharge model.

Outputs contain `[z,y,x]` moisture, pressure and total-head snapshots, a
material/active mask, hourly interval-average boundary flows, UTC times, a
signed water ledger and local nonlinear residual diagnostics. JSON uses
`null` for inactive cells. Backward Euler substeps are at most 600 s; hourly
output does not limit numerical stepping to one hour.

### Sensors and where their information enters

The [sensor guide](sensors.md) supplies complete CSV schemas and equations.
The committed fixture is synthetic, sampled from this same 72-hour run with
deterministic offsets. It demonstrates ingestion and held-out bookkeeping;
it does not supply independent physical truth.

| Example sensor | Forward quantity | Intended future role |
| --- | --- | --- |
| SM_1 and held-out SM_2 | Water-content average over declared support | State/parameter constraints and independent-site assessment |
| PSI_1 | Gauge pore pressure `rho g psi` | Slip-depth pressure/suction constraints; no strength or failure calculation |
| HEAD_1 | Uniform-conductance screen average in one saturated aquifer | Head/well-level constraint with a declared screen and datum |
| Q_OUT | Exact observation-window average of identified outlet flux | Integrated flow constraint once a suitable river operator exists |
| P_GAUGE | Precipitation rate/accumulation over its window | Uncertain forcing evaluation; not an independent state update |

The ingestion boundary is `hydrostate.observations`: metadata and raw records
are checked, units normalized, forward quantities sampled and residuals grouped
by kind, units and held-out role. `--as-of` filters by arrival time for replay
experiments. A future analysis update belongs **after** these operators and
**before** forecast propagation, with constitutive reconciliation and a
separate analysis water increment. No state update is implemented here.

To generate another clearly labelled synthetic fixture from a different run:

```bash
hydrostate sensors --model runs/valley3d/results.json \
  --make-demo runs/sensors/generated --output runs/sensors/generated_report.json
```

### Resolution and timing

The valley's local/global conservation checks pass. A 24-hour probe includes
the initial transient and the beginning of the first rainfall pulse:

```bash
python scripts/check_valley_refinement.py --hours 24 \
  --output runs/valley3d/refinement.json
```

| Run | Active cells | Net east-boundary outflow (m³) | Storage change (m³) |
| --- | ---: | ---: | ---: |
| Reference, maximum step 600 s | 724 | 4.36295 | 19.63705 |
| Same grid, maximum step 300 s | 724 | 4.36327 | 19.63673 |
| Double resolution in x/y/z, 600 s | 5,936 | 5.61762 | 18.38238 |

Halving the time step changes this cumulative outflow by **0.0073%**.
Doubling grid resolution changes it by **28.8%**, while increasing active
domain volume from 36,200 to 37,100 m³ and initial stored water by 245 m³.
The coarse-footprint surface elevation RMSE is 0.183 m. These are combined
geometry/material/solution sensitivities, not a clean discretization-order
estimate. Top-cell means also change their physical sampling elevation.
The [compact measured evidence](technical/valley_refinement.json) preserves
those differences and comparison support. All three budgets pass and the
boundary remains unponded.

The physical valley solution is **not established as mesh-converged**.
Stair geometry, the lens,
water-table position, surface pressure and boundary discharge must be checked
under horizontal/vertical refinement before using those outputs quantitatively.
The measured 24-hour runs took approximately 1.3 s, 2.6 s and 65 s respectively
on this macOS ARM environment, excluding imports and JSON writing. The current
sparse direct solver's runtime growth motivates profiling and a scalable
iterative/preconditioned solver; it does not identify a bottleneck without a
profile. Peak resident memory was not measured. These are small-domain
benchmarks, not watershed/Rainier/regional forecast-capacity estimates.
Record end-to-end
ingestion, ensemble, assimilation, restart and output costs for that claim.

## 3. Hillslope rain-on-snow and a postfire perturbation

```bash
python -m puget_modflow.toy_watershed --output runs/toy_watershed/results.json
python scripts/plot_toy_watershed.py --input runs/toy_watershed/results.json \
  --output runs/toy_watershed
python scripts/check_toy_refinement.py --help
```

This independent 2D reference adds snow/interception stores, stressed root
uptake, ponding and a lumped channel to mixed Richards flow. The scenarios
exercise process accounting and pressure exports. Prescribed degree-day melt
is not DHSVM's energy balance. The postfire parameter changes are hypothetical
sensitivity perturbations, not a universal fire-response law.

[Measured refinement](toy_watershed.md) found substantial vertical-grid
sensitivity in runoff and the outlet hydrograph. Those results remain visible;
neither a tight water budget nor the new 3D linear verification resolves them.

## 4. Next experiments and the evidence needed

| Experiment | Scientific question | Acceptance evidence before interpretation |
| --- | --- | --- |
| Layered column with moving water table | Can upward capillary supply and drainage coexist without duplicate storage? | Independent variably saturated reference; front/pressure and mesh/time convergence; moving-interface budgets |
| Refined 3D valley and saturated riparian zone | How do a finite lens and lateral gradients redistribute pressure? | Resolved geometry/layers; pressure and discharge convergence; below-/above-channel-bed seepage conditions |
| Wet valley with losing/gaining stream | Does routing respond consistently to groundwater feedback? | Surface/channel storage owner, flow reversal, dry bed and flood-stage tests; reach-wise exchange closure |
| Rooted dry-to-wet valley | How does water availability alter uptake and oxygen stress? | Separate transpiration/evaporation demands; root/sap-flow/profile data; compensation and redistribution benchmarks |
| Delayed-data replay | Do head/moisture/flow observations improve forecasts away from training sites? | Replay equivalence; analysis increments in ledger; held-out sites/events; forecast skill and uncertainty calibration |
| DHSVM–Richards–MODFLOW pilot | Is a nonoverlapping split accurate and fast enough for hourly cycles? | Conservatively mapped volumes, converged bidirectional coupling, complete restart and measured full-cycle latency |

Hazards use the hydrologic fields as antecedent states. Landslide assessment
adds stability/strength and slip-surface geometry; earthquake liquefaction and
lateral spreading add shaking and deformation/triggering physics; floods add
channel/floodplain hydraulics; wildfire adds fuel/weather/fire behavior.
Hourly antecedent states do not resolve seconds-long shaking or short
postfire rainfall bursts. See [AGENTS.md](../AGENTS.md) for the application
contract and [feasibility](coupling_feasibility.md) for scale-specific limits.
