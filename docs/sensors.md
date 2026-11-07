# Sensor ingestion and validation

`hydrostate sensors` reads portable CSV metadata and measurements, applies
explicit forward operators to a precomputed Cartesian 3-D model run, and writes
a JSON comparison report. It performs **no assimilation, calibration or model
state updates**. The committed examples are synthetic samples from the same
valley simulation, with deterministic small offsets and one missing record.
Their residuals test the pipeline; they do not establish field accuracy.

## Run the example

After generating `runs/valley3d/results.json`, run:

```bash
python scripts/validate_sensors.py \
  --model runs/valley3d/results.json \
  --sensors examples/sensors/sensors.csv \
  --observations examples/sensors/observations.csv \
  --output runs/sensors/report.json
```

The equivalent installed command is `hydrostate sensors` with the same flags.
To deliberately regenerate the synthetic CSVs from a new simulation:

```bash
python scripts/validate_sensors.py \
  --model runs/valley3d/results.json \
  --make-demo examples/sensors \
  --output runs/sensors/report.json
```

`--make-demo DIR` writes/replaces `DIR/sensors.csv` and `DIR/observations.csv`.
Normal comparison does not edit the input files. `--as-of` accepts an explicit
timezone-aware timestamp and excludes measurements arriving after that time.
Malformed metadata, ambiguous duplicates, unsupported operators and reference
mismatches terminate with exit status 2; inspect stderr and the exit status.
The JSON result contains input paths and SHA-256 hashes for provenance.

## CSV contract

Use UTF-8 with a header row. Empty optional fields remain empty. Coordinate
columns are always Cartesian metres, with `z_m` positive upward; measurement
units are specified separately. Source coordinates and measured elevations
must be transformed to the model's CRS/datum before ingestion. This adapter
checks exact reference labels and performs no coordinate or datum transform.

`sensors.csv` has one row per instrument:

| Columns | Meaning |
| --- | --- |
| `sensor_id`, `kind` | Unique stable ID and one supported quantity below |
| `x_m`, `y_m` | Model-coordinate location, including for precipitation/outlet records |
| `z_m` or `depth_m` | Exactly one for nonscreen state sensors; depth is vertical metres below the model surface in the containing horizontal cell |
| `support_x_m`, `support_y_m`, `support_z_m` | Full widths of an axis-aligned measurement support; empty/zero means point support in that direction |
| `screen_top_m`, `screen_bottom_m` | Optional well-screen elevations in the same vertical datum, top greater than bottom; replaces `z_m`/`depth_m`/vertical support |
| `aquifer_id`, `screen_assumption` | Required for screens; ID must match the model's declared aquifer, and assumption must equal `single_aquifer_uniform_conductance` |
| `outlet_id` | Required for discharge; must match `meta.outlet.id` rather than guessing a gauge from the nearest grid cell |
| `units`, `uncertainty` | Measurement units and positive 1-sigma uncertainty in those units |
| `horizontal_crs`, `vertical_datum` | Explicit references matching model metadata |
| `role` | `train` or `validation`; reports keep them separate |
| `provenance` | Source and processing history; demo rows explicitly identify synthetic origin |

Screens are supported for `hydraulic_head` and `water_level_depth`. They assume
one hydraulically connected aquifer, no pumping/wellbore-flow effects, and
uniform conductance per unit screen length. The code checks the declared
aquifer and nonnegative pressure in contributing model cells; a label does not
prove geological connectivity or that an observed screen is fully submerged.
Use short-screen/nested piezometers or a richer well observation model when
vertical gradients, confining layers, partial submergence or well flow matter.
In particular, nonnegative cell-centre pressure cannot rule out the water table
cutting the upper part of a screen within a cell. Resolve that crossing with
finer layers or a suitable subcell/well model before interpreting such a screen.

USGS monitoring guidance emphasizes surveyed horizontal/vertical references,
stable measuring points and hydraulic communication with the aquifer.
Convert a tape measurement below the top of casing to the stated land-surface
reference before using `water_level_depth` here.
[USGS monitoring guidance](https://pubs.usgs.gov/circ/2001/circ1217/html/essential.html),
[USGS measuring-point procedures](https://pubs.usgs.gov/tm/1a1/pdf/tm1-a1.pdf).

`observations.csv` has one row per sensor and valid time:

| Columns | Meaning |
| --- | --- |
| `sensor_id` | Existing metadata ID |
| `valid_time_utc` | Instant of a state measurement or end of a rate/accumulation interval |
| `arrival_time_utc` | Actual availability time, at or after valid time |
| `interval_start_utc` | Required for precipitation/discharge; empty for instantaneous states |
| `value` | Numeric observation; empty/NA/null remains missing, never zero-filled |
| `units` | Optional override of the sensor's measurement units |
| `uncertainty` | Optional 1-sigma override in the row's units; otherwise metadata uncertainty retains its metadata units |
| `quality_flag` | `good`, `provisional`, `suspect`, `bad` or `missing` |

All timestamps require an explicit timezone. UTC `Z` is recommended; explicit
offsets are normalized to UTC. Naive local dates are rejected, including dates
that would be ambiguous at a daylight-saving transition. Duplicate sensor/
valid-time pairs are rejected after normalization; resolve revisions upstream.
Only `good` observations enter default metrics. Missing values, excluded QC and
late arrivals are listed with reasons. Temporal averaging of state observations
is not implemented and must not be represented as an instantaneous sample.

## Quantities and units

| Kind | Model quantity | Canonical output | Accepted measurement units |
| --- | --- | --- | --- |
| `soil_moisture` | Volumetric `theta` over the stated support | `m3/m3` | `m3/m3`, `%` and common cubic-metre aliases |
| `pressure_head` | `pressure_head_m`, signed | `m` | `m`, `cm`, `mm`, `ft` |
| `pore_pressure` | `1000 × 9.80665 × pressure_head_m` | `Pa` | `Pa`, `kPa` |
| `hydraulic_head` | `total_head_m` at point/screen | `m` | `m`, `cm`, `mm`, `ft` |
| `water_level_depth` | Local surface elevation minus hydraulic head | `m` | `m`, `cm`, `mm`, `ft` |
| `streamflow` | Identified outlet's interval-mean `outlet_m3_s` | `m3/s` | `m3/s`, `L/s`, `cfs` |
| `precipitation` | Uniform forcing's interval mean `precip_mm_h` | `mm/h` | `mm/h`, `m/s`, `in/h`, or `mm`/`in` accumulated over the explicit interval |

Negative pressure denotes matric suction relative to atmospheric pressure.
The pressure conversion assumes freshwater density 1000 kg/m³ and standard
gravity 9.80665 m/s²; it is not a variable-density or absolute-pressure operator.
Negative `water_level_depth` is an artesian level above the surface and is not
clipped. A confined well's potentiometric level is not an unconfined water
table. No water-table crossing or perched-zone identification is inferred from
a head-well measurement.

Moisture observations must be volumetric rather than gravimetric; percentages
are divided by 100. `ft` uses the international foot (0.3048 m) and `cfs` the
corresponding cubic foot. Historic survey-foot datasets require explicit
preprocessing; datum conversion is separate from unit conversion.
[NIST unit definitions](https://www.nist.gov/pml/us-surveyfoot/revised-unit-conversion-factors).

## Forward-operator mathematics

For an interior point, trilinear weights in Cartesian cell-centre coordinates
give `H[f] = Σ w_kji f_kji`, with weights summing to one. In the half-cell
between the outermost centre and the grid edge, the point operator uses that
boundary cell value rather than extrapolating. Rectangular support uses
cell-overlap fractions for each nonzero support dimension, treating model cell
values as piecewise constant there. This is a representation choice, not an
instrument-response calibration. A real probe footprint may require a more
specific kernel and representativeness uncertainty.

The screen operator uses `H[h] = Σ L_k h_k / Σ L_k`, where `L_k` is the length
of the screen intersecting a vertical model cell. Horizontal location uses the
same interpolation rule. It is the constant-conductance-per-unit-length limit
of a screen-weighted head observation. Both partially inactive support and
partially unsaturated contributing screen cells are rejected; support is not
silently moved or renormalized onto wet cells. Bounds, masked cells and missing
model values are errors. Interior state times between stored snapshots use
linear time interpolation; no extrapolation or model rerun is performed.

Rate records at `t_i` describe the preceding model interval `(t_(i-1), t_i]`.
For observation interval `(a,b]`, the mean is

```text
H[q] = sum_i(q_i * overlap((a,b], (t_(i-1),t_i])) / (b-a).
```

The initial record is not treated as an earlier interval. Partial-hour and
multi-hour intervals are integrated by overlap; rates are not linearly
interpolated, shifted by half a step or sampled at their endpoint. Rain-gauge
comparisons require `meta.precipitation_spatial_support = "uniform"`; a domain
mean alone is insufficient to predict a point gauge under nonuniform forcing.
Precipitation metrics are separate **forcing comparisons**, not evidence of
state assimilation. Shared gauge/forcing errors still need explicit treatment
when real data are used.

The current valley outlet is an identified fixed-head regional boundary flux
surrogate. Comparing its synthetic `Q_OUT` record does not validate river
routing, reach gain/loss, stage, discharge ratings or flood hydraulics.

## Metrics and independent validation

Reports group by **quantity, canonical unit and role**, with per-sensor results
and matched/excluded records. For residual `r = model − observation`, they
include bias, RMSE, MAE, inverse-observation-variance weighted RMSE, mean
`(r/sigma_obs)^2`, and the fraction inside `±1.96 sigma_obs`. These are
descriptive statistics. Observation-only sigma does not create a forecast
confidence interval or a calibrated innovation chi-square: that requires model
uncertainty, bias and error-covariance assumptions. Correlated repeated sensors
and forcing products are not automatically independent samples.

`role=validation` keeps a site out of the `train` summary; this module does not
fit either set. The investigator must establish that validation sites/events
were withheld from parameter selection, forcing correction and tuning.
Field validation also needs calibrated instruments, surveyed references,
representativeness estimates and multiple independent events. A same-model
synthetic truth with small residuals is only a software self-consistency check.

Analytical tests in `tests/test_observations.py` cover linear spatial/time
fields, finite support, screen weights, pressure-versus-head identity,
piecewise-constant interval integration, units/uncertainty, UTC normalization,
duplicates, QC/missing values, references/bounds/inactive cells, holdouts,
known residuals and a complete CSV/JSON CLI round trip. Run them with:

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -p test_observations.py -v
```
