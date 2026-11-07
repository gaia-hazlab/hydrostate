# Synthetic sensor example

These CSVs are **synthetic**, sampled from the project's 72-hour Cartesian 3-D
valley run with small deterministic offsets. They are not downloaded instrument
data, independent truth, or evidence of field validation. `SM_1` has the `train`
role; the other sites have `validation` roles. One empty `SM_1` measurement is
explicitly flagged missing and must not become zero.

The example includes volume-supported soil moisture in percent, a signed pore
pressure in kPa, a screened hydraulic head in metres, an outlet boundary-flux
surrogate in L/s, and a uniform precipitation forcing comparison in mm/h.
Screens refer to one declared synthetic aquifer; arbitrary local elevations
must not be relabelled as NAVD88 or another real datum.

Regenerate deliberately after running the valley model:

```bash
python scripts/validate_sensors.py --model runs/valley3d/results.json \
  --make-demo examples/sensors --output runs/sensors/report.json
```

To use the committed files without replacing them, use `--sensors sensors.csv`
and `--observations observations.csv` with their full paths. See
[`docs/sensors.md`](../../docs/sensors.md) for the schema, operators, uncertainty
meaning, supported units and limits. Reports retain input SHA-256 hashes.
