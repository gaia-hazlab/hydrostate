# Legacy Puget Sound pipeline

This guide preserves the original monthly pipeline. The main project is now Hydrostate; see the [project README](../README.md) for the 3D experiments and verification.


This repository scaffolds a reproducible workflow for predicting groundwater
levels in the Puget Sound lowland from precipitation forcing and the best
available public hydrogeologic information for western Washington.

The default scientific strategy is conservative:

1. Use the USGS southeast Puget Sound hydrogeologic framework and MODFLOW-NWT
   model as the reference model where the study area overlaps King and Pierce
   Counties.
2. Use PRISM gridded precipitation as the precipitation forcing.
3. Prefer the published USGS Soil-Water Balance recharge model for recharge
   estimation; use the included precipitation-to-recharge transform only for
   scenario screening or software smoke tests.
4. Calibrate and validate against USGS NWIS and Washington Ecology groundwater
   level observations before making management-grade predictions.

The included "template" MODFLOW 6 builder is not a calibrated Puget Sound
groundwater model. It exists to verify that the software pipeline can write and
run MODFLOW inputs. For scientific predictions, ingest the USGS reference model
or derived calibrated package inputs before running scenarios.

## Repository Layout

```text
config/puget_sound.toml       Default pipeline configuration
docs/data_sources.md          Source notes and citations
src/puget_modflow/            Pipeline package and CLI
tests/                        Dependency-light unit tests
Makefile                      Common local commands
```

## Quick Start

Create the Pixi environment with the geospatial, MODFLOW, calibration, and test
dependencies:

```bash
pixi install
```

Run the lightweight tests:

```bash
pixi run test
```

Check that FloPy can write a MODFLOW 6 workspace and that the Pixi-provided
`mf6` executable can run it:

```bash
pixi run build-template-smoke
pixi run run-template
```

If you prefer a plain virtual environment instead, install the Python package
directly:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[modflow,geo,calibration,dev]"
```

Install MODFLOW executables separately from USGS and make sure `mf6` is on
`PATH`. If you plan to run the published USGS reference model directly, also
install the MODFLOW-NWT executable used by that model.

Initialize or refresh the default configuration:

```bash
pixi run python -m puget_modflow init-config --output config/puget_sound.toml
```

Write a concrete execution plan:

```bash
pixi run plan
```

Fetch ScienceBase metadata for the reference datasets:

```bash
pixi run fetch-reference
```

Add `--download-files` when you are ready to download the attached model and
support archives; those files can be large.

Download monthly PRISM precipitation grids for the calibration period:

```bash
pixi run python -m puget_modflow fetch-prism \
  --config config/puget_sound.toml \
  --start 2005-01 \
  --end 2015-12 \
  --output data/raw/prism/monthly_ppt
```

Aggregate the PRISM grids to a domain-mean monthly precipitation CSV. This step
requires the `geo` optional dependencies:

```bash
pixi run python -m puget_modflow aggregate-prism \
  --config config/puget_sound.toml \
  --input data/raw/prism/monthly_ppt \
  --output data/processed/monthly_precip.csv
```

If you have already aggregated PRISM or station precipitation to monthly basin
means, create a CSV with columns `date` and either `precip_in` or `precip_mm`,
then build a uniform recharge time series:

```bash
pixi run python -m puget_modflow recharge-from-csv \
  --config config/puget_sound.toml \
  --precip-csv data/processed/monthly_precip.csv \
  --output data/processed/recharge_timeseries.csv
```

Build a small MODFLOW 6 template model for a software smoke test:

```bash
pixi run build-template
```

Run a MODFLOW workspace:

```bash
pixi run python -m puget_modflow run --workspace runs/template_mf6 --exe mf6
```

## Prediction Workflow

For a scientific Puget Sound prediction run:

1. Download the USGS ScienceBase datasets listed in `docs/data_sources.md`.
2. Select the model domain:
   - `southeast_puget_sound` for the published 13-layer USGS active model area.
   - `regional_puget_lowland` for exploratory regional studies that need a new
     grid and a coarser hydrogeologic framework.
3. Build recharge using the USGS SWB model when possible. If only precipitation
   is available, calibrate recharge fractions against observed groundwater
   levels and baseflow before using the scenario outputs.
4. Run steady-state baseline, transient historical calibration, and prediction
   scenarios.
5. Compare simulated heads to USGS NWIS and Washington Ecology well data, and
   report prediction uncertainty.

## Tests

The core tests do not require FloPy or geospatial packages:

```bash
PYTHONPATH=src python3 -m unittest discover -s tests
```

## Vadose-Zone Research Prototype

The proposed coupled architecture and its process/validation requirements are
described in [the vadose design](vadose_design.md), with separate
[spatial parallelization and optimization notes](spatial_parallelism.md).
Project-wide physical and hazard requirements are in [AGENTS.md](../AGENTS.md).

An independent synthetic hillslope-and-channel demonstrator exercises mixed-form
Richards flow in two dimensions, layered soil, capillarity, lateral flow, root
uptake, surface storage and channel discharge. It uses hypothetical parameters
and simplified snow/interception processes. Its saturated and unsaturated
regions share one solver; **it is not yet a DHSVM–MODFLOW coupling**, a
calibrated watershed, or a hazard prediction model. Its pressure fields are
hydrologic inputs for future ecological/geotechnical analyses.

The outlet hydrograph is **not mesh-converged**: refining the vertical grid
materially changes infiltration, ponding and runoff. See the
[toy setup, measured checks and limitations](toy_watershed.md) before
interpreting any output quantitatively.

Run and export the toy with the project environment:

```bash
pixi run toy-watershed
pixi run plot-toy
```

Alternatively, install the optional numerical dependencies with
`python -m pip install -e ".[toy]"`, then run:

```bash
python -m puget_modflow.toy_watershed --output runs/toy_watershed/results.json
python scripts/plot_toy_watershed.py --input runs/toy_watershed/results.json --output runs/toy_watershed
```

The results include hourly states and diagnostics in JSON, a CSV time series,
and PNG/PDF summary figures under `runs/toy_watershed/`. The numerical tests run
when NumPy/SciPy are installed and otherwise skip cleanly, preserving the
dependency-light pipeline test suite. A passing water budget is a numerical
check; independent reference comparisons and mesh/time refinement remain
separate validation requirements.
