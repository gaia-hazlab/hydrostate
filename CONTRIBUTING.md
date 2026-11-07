# Contributing to Hydrostate

Read [AGENTS.md](AGENTS.md) before changing physical equations, coupling,
observation operators or numerical algorithms. Describe the scientific problem,
applicable regime, governing equations and parameter provenance alongside the
implementation. Small, reviewable changes with reproducible experiments are
preferred.

## Local development

```bash
python -m pip install -e ".[toy,dev]"
python -m unittest discover -s tests
python -m hydrostate.verification --output runs/verification/report.json
```

The numerical test suite requires NumPy/SciPy; absence of those dependencies
must not masquerade as a passing numerical verification run. CI installs them.
Raw data, environments and generated run directories are ignored. Commit small
synthetic fixtures, configurations and compact reference metrics instead of
large generated fields or private sensor records.

## Physical and mathematical changes

Explain units, signs, geometry, storage ownership and time support. Test a
known solution or independently derived property in addition to water closure.
Use analytical/manufactured cases to measure errors and convergence orders.
For nonlinear flow, retain independent infiltration/capillary benchmarks and
mesh/time sensitivity; a linear manufactured case does not verify every
nonlinear regime.

Run relevant unit tests and experiments. Record solver tolerances, software
version/commit, hardware and measured runtime when making performance claims.
Optimization comparisons must retain the same equations, forcing, grid and
acceptance criteria; report approximation changes separately.

## Observations and scientific claims

Keep sensor support, units, datum/CRS, quality, uncertainty, valid time and
arrival time explicit. Synthetic observations test interfaces and observation
operators; they do not establish field validity. Withhold sensors/events from
calibration or assimilation before claiming prediction skill.

Distinguish implemented features, mathematically checked regimes, unresolved
numerical sensitivities and field validation. Hydrologic conditioning alone
does not establish liquefaction, landslide, flood inundation or wildfire skill.
Report new limitations directly in the relevant technical report.
