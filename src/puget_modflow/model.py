from __future__ import annotations

from pathlib import Path
import csv
import json
import shutil
import subprocess

from puget_modflow.config import PipelineConfig
from puget_modflow.recharge import read_recharge_timeseries


def write_run_manifest(
    config: PipelineConfig,
    output_dir: str | Path,
    *,
    extra: dict[str, object] | None = None,
) -> Path:
    workspace = Path(output_dir)
    workspace.mkdir(parents=True, exist_ok=True)
    payload = config.to_dict()
    if extra:
        payload["extra"] = extra
    manifest = workspace / "run_manifest.json"
    manifest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return manifest


def build_mf6_template(
    config: PipelineConfig,
    output_dir: str | Path,
    *,
    recharge_csv: str | Path | None = None,
    exe_name: str = "mf6",
) -> Path:
    """Write a small MODFLOW 6 model for pipeline smoke testing.

    This is not the calibrated Puget Sound model. It confirms that FloPy can
    write packages, recharge forcing can be mapped to stress periods, and an
    installed `mf6` executable can run the workspace.
    """
    flopy = _require_flopy()
    import numpy as np

    workspace = Path(output_dir)
    workspace.mkdir(parents=True, exist_ok=True)

    recharge = _load_recharge_by_period(recharge_csv)
    nper = max(len(recharge), 1)
    perioddata = [(30.4375, 1, 1.0) for _ in range(nper)]

    sim = flopy.mf6.MFSimulation(
        sim_name=config.project.name,
        version="mf6",
        exe_name=exe_name,
        sim_ws=str(workspace),
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units=config.model.time_units.upper(),
        nper=nper,
        perioddata=perioddata,
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        outer_dvclose=1e-4,
        inner_dvclose=1e-4,
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname="puget_template",
        save_flows=True,
    )
    botm = np.linspace(
        config.model.top - (config.model.top - config.model.bottom) / config.model.nlay,
        config.model.bottom,
        config.model.nlay,
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=config.model.length_units.upper(),
        nlay=config.model.nlay,
        nrow=config.model.nrow,
        ncol=config.model.ncol,
        delr=config.model.delr,
        delc=config.model.delc,
        top=config.model.top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=config.model.top - 25.0)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=[1] + [0] * (config.model.nlay - 1),
        k=25.0,
        k33=2.5,
    )
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=[1] + [0] * (config.model.nlay - 1),
        ss=1e-5,
        sy=0.12,
        steady_state={0: False},
        transient={0: True},
    )

    first_layer_bottom = float(botm[0])
    left_head = (config.model.top + first_layer_bottom) / 2.0
    right_head = config.model.top - (config.model.top - first_layer_bottom) * 0.25
    chd_data = {}
    for period in range(nper):
        cells = []
        for row in range(config.model.nrow):
            cells.append(((0, row, 0), left_head))
            cells.append(((0, row, config.model.ncol - 1), right_head))
        chd_data[period] = cells
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_data, pname="CHD")

    rch = {period: recharge.get(period, 1e-4) for period in range(nper)}
    flopy.mf6.ModflowGwfrcha(gwf, recharge=rch, pname="RCHA")

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord="puget_template.cbc",
        head_filerecord="puget_template.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    sim.write_simulation()
    write_run_manifest(
        config,
        workspace,
        extra={
            "template_warning": "This is a software smoke-test model, not a calibrated Puget Sound prediction model.",
            "recharge_csv": str(recharge_csv) if recharge_csv else None,
        },
    )
    return workspace


def run_modflow_workspace(
    workspace: str | Path,
    *,
    exe_name: str = "mf6",
    timeout: int | None = None,
) -> subprocess.CompletedProcess[str]:
    executable = shutil.which(exe_name)
    if executable is None:
        raise FileNotFoundError(f"MODFLOW executable not found on PATH: {exe_name}")
    return subprocess.run(
        [executable],
        cwd=Path(workspace),
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )


def _load_recharge_by_period(path: str | Path | None) -> dict[int, float]:
    if path is None:
        return {0: 1e-4}
    records = read_recharge_timeseries(path)
    return {record.stress_period: record.recharge_ft_per_day for record in records}


def _require_flopy():
    try:
        import flopy
    except ImportError as exc:
        raise RuntimeError(
            "FloPy is required to build MODFLOW workspaces. Install with "
            '`python -m pip install -e ".[modflow]"`.'
        ) from exc
    return flopy


def write_head_observation_template(path: str | Path) -> Path:
    """Write the CSV schema expected for calibration observation targets."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "site_no",
        "site_name",
        "date",
        "x",
        "y",
        "layer",
        "observed_head_ft",
        "screen_top_ft",
        "screen_bottom_ft",
        "source",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
    return output
