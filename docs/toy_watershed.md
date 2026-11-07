# Synthetic watershed demonstrator

This is a working two-dimensional mixed-Richards numerical prototype, with a
small idealized headwater contributing area and lumped channel. It implements
a subset of the [target vadose/ecological/geotechnical design](vadose_design.md).
It does not invoke DHSVM or MODFLOW, assimilate observations, or predict a
hazard. All geometry, soil, vegetation and forcing parameters are hypothetical.

**The outlet hydrograph is not mesh-converged.** Water-budget and solver checks
pass, but refinement materially changes runoff. Use this example to inspect
states, conservation, numerical sensitivity and software interfaces, not as
quantitative watershed or geotechnical prediction.

## Implemented model

- A 120 m downslope section, 6 m normal depth and 100 m out-of-plane width,
  contributing approximately 1.2 ha in horizontal projection. Slope is 3%.
  The 24 × 24 finite-volume grid has 5 m downslope spacing and 0.25 m layers.
- The coordinates form an orthogonally rotated rectangle. Actual elevation is
  `z = z0 - s sin(beta) - d cos(beta)`. Atmospheric vertical precipitation is
  converted from horizontal projected area to slope-area flux using `cos(beta)`.
- Van Genuchten–Mualem retention/conductivity, a low-conductivity lens between
  1.5 and 2.25 m depth over part of the slope, and weathered substrate below
  2.5 m. Capillarity, lateral/vertical redistribution and saturated pressure
  occur in one continuum; shallow groundwater is not duplicated.
- Cell-integrated mixed storage, backward Euler and modified Picard iterations.
  The declared toy storage law is `theta(psi) + Ss max(psi,0)`, with small
  saturated elastic storage. Acceptance checks both volume-rate residuals and
  pressure increments. Rejected steps restore the beginning state.
- A dry/ponded top-boundary switch, implicit local surface detention and a
  channel reservoir. Wet bank exchange is bidirectional and restricted to the
  actual bank above the channel floor. Exposed bank faces allow seepage;
  the below-floor toe wall and the divide/base are impermeable. Horizontal
  streambed exchange and spatial overland hydraulics are omitted.
- Synthetic snowfall/rain, a conservative snow store with temperature-index
  melt, canopy interception/evaporation and root uptake limited by dry and
  near-saturated conditions. The wet-end factor is an empirical uptake limit,
  not an oxygen-transport calculation. `root_stress` is an uptake multiplier:
  one means unstressed and zero means shut down.
- Hourly output with adaptive internal steps up to 600 seconds. Records are
  end-of-hour states; rate variables represent means over the preceding hour.
  The artificial start is 2026-01-01T00:00:00Z, not an observed winter event.

The basic references are [van Genuchten, 1980](https://doi.org/10.2136/sssaj1980.03615995004400050002x)
and [Celia et al., 1990](https://doi.org/10.1029/WR026i007p01483). This toy is
inspired by their constitutive/mixed-form methods; it does not reproduce an
independent published validation benchmark.

## Forcing and outputs

The baseline begins with 25 mm SWE. It adds 0.8 mm/h precipitation over hours
6–18 under freezing conditions, then warm rain of 7 mm/h over hours 30–42,
12 mm/h over 60–66, and 3 mm/h over 85–90. The warmer phase melts stored snow.
These windows are half-open; hourly output rates are timestamped at their ends.
Potential ET has a prescribed diurnal cycle, allocated first to interception
evaporation and then to roots.

The default run spans 120 hours and exports:

- `results.json`: hourly fields and scalar states, grid/material labels,
  parameter/forcing assumptions, water ledger and solver diagnostics.
- `summary.png` and `summary.pdf`: initial/wettest moisture sections, zero
  pressure contours, forcing/outlet trace, root moisture and shallow pressure.
- `timeseries.csv` and `summary.json`: scalar data and run summary.
- `refinement.json`: half-time-step and doubled-grid comparisons.

The pressure sensor is at downslope coordinate 62.5 m and normal depth 1.125 m.
Negative reported pore pressure is matric suction. The water-table diagnostic
tracks a base-connected zero-pressure crossing and excludes isolated perches.
A zero-pressure contour alone does not identify one regional water table.

## Measured numerical checks

The baseline completed 120 simulated hours in approximately 5.8 seconds in the
project environment. This is a tiny serial synthetic problem; the figure
excludes plotting/startup and establishes no operational runtime for DHSVM,
MODFLOW, real terrain, data assimilation or ensembles. The 48 × 48 comparison
took about 18 seconds in a separate run.

| Baseline quantity | Result |
| --- | ---: |
| Precipitation | 2,166.2254 m³ |
| ET, canopy plus roots | 81.7558 m³ |
| Channel outlet volume | 201.0616 m³ |
| Total storage increase, including snow/canopy/channel | 1,883.4080 m³ |
| Signed total budget residual | −1.38 × 10⁻⁷ m³ |
| Equivalent residual depth | −1.15 × 10⁻⁸ mm |
| Maximum accepted local mass residual | 1.98 × 10⁻⁹ m³/s |
| Maximum accepted pressure increment | 9.96 × 10⁻⁶ m |

The budget is `final storage - initial storage - precipitation + ET + outlet`.
Internal melt, infiltration, runoff and bank exchange cancel; initial snow is
included in initial storage. Normalized error divides the absolute residual by
`max(total precipitation, total ET + outlet)`, excluding initial storage.

The baseline root-zone water content changes from about 0.224 to a maximum
0.389 m³/m³, and the shallow sensor's pressure changes from approximately
−23.6 toward −3.6 kPa. These illustrate wetting and loss of suction in the
chosen discretization; they are not measured site conditions.

| Sensitivity comparison | Maximum root-moisture difference | Storm peak difference | Total outlet-volume difference |
| --- | ---: | ---: | ---: |
| Maximum step 600 → 300 s | 0.000213 m³/m³ | +3.0% | +1.2% |
| Grid 24 × 24 → 48 × 48 | 0.00531 m³/m³ | −86.4% | −35.5% |

Storm peaks here are evaluated after hour 24 to exclude initial reservoir
recession. The finer-grid storm peak occurs 29 hours later; it is a different
dominant pulse. Final pressure-head RMSE on common reference locations is
approximately 0.0017 m for the time comparison and 0.136 m for the grid
comparison. Complete metrics and normalization are in `refinement.json`.

Both refined runs close their budgets. The large spatial sensitivity remains
material. Mesh volumes, lens bounds and initial bank area were checked to be
consistent. Unresolved wetting/ponding and nonlinear face-conductivity treatment
are candidate causes, not an established diagnosis. Directional refinement
and independent infiltration/ponding benchmarks are required before using the
hydrograph quantitatively. Do not tune conductivity averaging just to force
agreement between these two runs.

Directional experiments in `resolution_diagnosis.json` identify vertical
resolution as the dominant sensitivity: doubling only vertical layers reduces
the storm peak by approximately 89%, whereas doubling only horizontal cells
changes it by approximately 0.12%. This narrows the next benchmark to vertical
infiltration/front and ponding treatment; it does not independently validate
the proposed mechanism or the finer solution.

The eight new numerical tests cover constitutive bounds/derivatives,
hydrostatic equilibrium on the rotated grid, closed redistribution, wet-bank
losing exchange and its ledger, below-floor no-flow, storm storage closure,
time-step sensitivity and excluding perches from the regional water-table
diagnostic. The complete project suite has 23 tests with numerical dependencies
installed. These checks are distinct from independent physical validation.

## Scope and next engineering work

Bare-soil/pond evaporation, compensated uptake and plant hydraulic
redistribution are omitted. Also absent are snow energy balance, glaciers,
heat/vapor/ice, hysteresis, macropores/fractures, deformation, geotechnical
triggering, fuel moisture/fire spread and flood/debris-flow hydraulics. The
target design specifies when and how to add these processes.

The channel begins with stored water; the first hours include unspun-up
reservoir recession. Supported scenarios keep the channel wet. Complete drying
requires supply-limited dry-bank handling and currently causes an explicit
failure. The `postfire` option is one assumed sensitivity branch, not a
universal prediction of fire effects; it reduces topsoil conductivity,
interception and vegetation demand. An `intense_storm` option doubles the
later storm intensity. All three supported scenarios completed with closed
budgets.

The current kernel precomputes face connectivity but assembles sparse matrices
and work arrays during iteration. After independent benchmarks, useful
engineering experiments include CSR pattern reuse, persistent work buffers,
compiled constitutive/face assembly and spatial tile ownership, followed by
robust global preconditioning. The [spatial design](spatial_parallelism.md)
distinguishes vertical-block preconditioning of full lateral flow from the
physical approximation of independent columns. MPI execution is not implemented
in this toy.

## Reproduce

```bash
pixi run toy-watershed
pixi run plot-toy
```

For a direct project-environment invocation without environment resolution:

```bash
PYTHONPATH=src .pixi/envs/default/bin/python -m puget_modflow.toy_watershed --output runs/toy_watershed/results.json
.pixi/envs/default/bin/python scripts/plot_toy_watershed.py --input runs/toy_watershed/results.json --output runs/toy_watershed
PYTHONPATH=src .pixi/envs/default/bin/python scripts/check_toy_refinement.py
PYTHONPATH=src .pixi/envs/default/bin/python -m unittest discover -s tests
```

The optional standalone numerical dependencies are installed by
`python -m pip install -e ".[toy]"`. The CLI supports `--scenario baseline`,
`--scenario intense_storm`, `--scenario postfire`, `--hours` and `--max-step-s`.
Numerical test classes skip only when NumPy/SciPy are absent; other import
errors remain visible failures.
