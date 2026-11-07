from __future__ import annotations

from pathlib import Path
import json

from puget_modflow.config import PipelineConfig
from puget_modflow.datasources import source_catalog


def build_pipeline_plan(config: PipelineConfig) -> dict[str, object]:
    workspace = config.project.workspace
    return {
        "project": config.project.name,
        "workspace": workspace,
        "domain": config.domain.__dict__,
        "strategy": [
            "Fetch authoritative USGS ScienceBase metadata and model-support files.",
            "Download precipitation forcing from PRISM for the calibration or scenario period.",
            "Aggregate precipitation to the model grid or active model area.",
            "Estimate recharge with the published SWB workflow when available; otherwise use calibrated precipitation fractions.",
            "Build or update MODFLOW packages.",
            "Run steady-state, transient historical, and prediction scenarios.",
            "Validate simulated heads against NWIS and Ecology observation wells.",
            "Export head maps, hydrographs, residuals, and uncertainty notes.",
        ],
        "commands": [
            {
                "name": "fetch_reference_metadata",
                "command": (
                    "python -m puget_modflow fetch-reference "
                    "--config config/puget_sound.toml --output data/raw/usgs_reference"
                ),
            },
            {
                "name": "fetch_prism_monthly",
                "command": (
                    "python -m puget_modflow fetch-prism "
                    "--config config/puget_sound.toml --start 2005-01 --end 2015-12 "
                    "--output data/raw/prism/monthly_ppt"
                ),
            },
            {
                "name": "aggregate_prism_to_monthly_precip",
                "command": (
                    "python -m puget_modflow aggregate-prism "
                    "--config config/puget_sound.toml "
                    "--input data/raw/prism/monthly_ppt "
                    "--output data/processed/monthly_precip.csv"
                ),
            },
            {
                "name": "build_recharge_from_aggregated_precip",
                "command": (
                    "python -m puget_modflow recharge-from-csv "
                    "--config config/puget_sound.toml "
                    "--precip-csv data/processed/monthly_precip.csv "
                    "--output data/processed/recharge_timeseries.csv"
                ),
            },
            {
                "name": "build_template_model",
                "command": (
                    "python -m puget_modflow build-template "
                    "--config config/puget_sound.toml "
                    "--recharge data/processed/recharge_timeseries.csv "
                    "--output runs/template_mf6"
                ),
            },
            {
                "name": "run_modflow",
                "command": "python -m puget_modflow run --workspace runs/template_mf6 --exe mf6",
            },
        ],
        "sources": source_catalog(config),
    }


def write_pipeline_plan(config: PipelineConfig, path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(build_pipeline_plan(config), indent=2), encoding="utf-8")
    return output
