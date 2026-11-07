"""Reproduce 3D valley time/mesh sensitivity with explicit staircase geometry.

Default 24 hours covers the initial transient and first rainfall pulse only.
Two grids do not establish convergence; the finer staircase also changes the
discretized physical domain. Storage changes are compared separately from the
different initial storage volumes. All three runs exercise the real kernel.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.interpolate import RegularGridInterpolator

from hydrostate.valley import ValleyConfig, run_valley


def summarize(result: dict) -> dict:
    grid = result["grid"]
    active = np.asarray(grid["active"], dtype=bool)
    psi = np.asarray(result["snapshots"][-1]["pressure_head_m"], dtype=float)
    theta = np.asarray(result["snapshots"][-1]["theta"], dtype=float)
    iz = active.sum(axis=0) - 1
    iy, ix = np.indices(iz.shape)
    volume = grid["dx_m"] * grid["dy_m"] * grid["dz_m"]
    budget = result["budget"]
    return {
        "config": result["meta"]["config"],
        "active_cells": int(active.sum()), "active_domain_volume_m3": float(active.sum() * volume),
        "mean_stair_surface_m": float(np.mean(grid["surface_elevation_m"])),
        "initial_storage_m3": budget["initial_storage_m3"],
        "final_storage_m3": budget["final_storage_m3"],
        "storage_change_m3": budget["final_storage_m3"] - budget["initial_storage_m3"],
        "precip_m3": budget["precip_m3"], "outlet_net_m3": budget["outlet_net_m3"],
        "final_domain_mean_theta": float(theta[active].mean()),
        "final_top_mean_theta": float(theta[iz, iy, ix].mean()),
        "final_top_mean_pressure_head_m": float(psi[iz, iy, ix].mean()),
        "signed_budget_residual_mm": budget["signed_residual_mm"],
        "budget_accepted": budget["accepted"], "diagnostics": result["diagnostics"],
    }


def compare(reference: dict, candidate: dict) -> dict:
    base, other = summarize(reference), summarize(candidate)
    coarse = reference["grid"]
    fine = candidate["grid"]
    z, y, x = np.meshgrid(coarse["z_m"], coarse["y_m"], coarse["x_m"], indexing="ij")
    query = np.column_stack((z.ravel(), y.ravel(), x.ravel()))
    candidate_pressure = np.asarray(candidate["snapshots"][-1]["pressure_head_m"], dtype=float)
    same_grid = all(np.array_equal(coarse[key], fine[key]) for key in ("x_m", "y_m", "z_m", "active"))
    if same_grid:
        # Exact centres need no interpolation. In particular, 0*NaN terms in
        # multilinear stencils must not discard valid cells beside the surface.
        interpolated = candidate_pressure
        support = "Direct arrays on identical coordinates and active mask; every active cell compared"
    else:
        interpolated = RegularGridInterpolator(
            (fine["z_m"], fine["y_m"], fine["x_m"]), candidate_pressure,
            bounds_error=False, fill_value=np.nan,
        )(query).reshape(z.shape)
        support = "Candidate trilinearly interpolated to reference cell centres; omit any stencil touching inactive cells; geometry changes remain in the comparison"
    base_pressure = np.asarray(reference["snapshots"][-1]["pressure_head_m"], dtype=float)
    common = np.isfinite(interpolated) & np.isfinite(base_pressure)
    pressure_difference = interpolated[common] - base_pressure[common]
    base_surface = np.asarray(coarse["surface_elevation_m"])
    fine_surface = np.asarray(fine["surface_elevation_m"])
    if fine_surface.shape == base_surface.shape:
        mean_fine_surface = fine_surface
    else:
        ny, nx = base_surface.shape
        mean_fine_surface = fine_surface.reshape(ny, 2, nx, 2).mean(axis=(1, 3))
    outlet_scale = abs(base["outlet_net_m3"])
    return {
        "candidate": other,
        "net_outlet_difference_m3": other["outlet_net_m3"] - base["outlet_net_m3"],
        "net_outlet_relative_difference": (other["outlet_net_m3"] - base["outlet_net_m3"]) / outlet_scale if outlet_scale else None,
        "storage_change_difference_m3": other["storage_change_m3"] - base["storage_change_m3"],
        "initial_storage_difference_m3": other["initial_storage_m3"] - base["initial_storage_m3"],
        "active_domain_volume_difference_m3": other["active_domain_volume_m3"] - base["active_domain_volume_m3"],
        "stair_surface_coarse_footprint_rmse_m": float(np.sqrt(np.mean((mean_fine_surface - base_surface) ** 2))),
        "top_mean_theta_difference": other["final_top_mean_theta"] - base["final_top_mean_theta"],
        "top_mean_pressure_difference_m": other["final_top_mean_pressure_head_m"] - base["final_top_mean_pressure_head_m"],
        "common_active_pressure_cells": int(common.sum()),
        "reference_active_cells": int(np.isfinite(base_pressure).sum()),
        "pressure_head_rmse_on_common_reference_centres_m": float(np.sqrt(np.mean(pressure_difference ** 2))),
        "pressure_support": support,
    }


def source_manifest() -> dict:
    root = Path(__file__).resolve().parents[1]
    files = ("src/hydrostate/flow.py", "src/hydrostate/valley.py", "scripts/check_valley_refinement.py")
    return {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in files}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hours", type=int, default=24)
    parser.add_argument("--output", type=Path, default=Path("runs/valley3d/refinement.json"))
    args = parser.parse_args(argv)
    c = ValleyConfig(hours=args.hours)
    print("Running reference 3D valley", flush=True)
    reference = run_valley(config=c, snapshot_interval_h=c.hours)
    print("Running half-time-step 3D valley", flush=True)
    time_refined = run_valley(config=replace(c, max_step_s=c.max_step_s / 2), snapshot_interval_h=c.hours)
    print("Running doubled x/y/z grid 3D valley", flush=True)
    mesh_refined = run_valley(config=replace(c, nx=2*c.nx, ny=2*c.ny, nz=2*c.nz), snapshot_interval_h=c.hours)
    result = {
        "interpretation": "Restricted synthetic nonponding initial-transient sensitivity; two grids do not establish convergence or field skill",
        "geometry_caveat": "Refining the staircase changes physical active volume, the surface and material representation; top-cell means also move in space",
        "software": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform(), "machine": platform.machine()},
        "source_sha256": source_manifest(),
        "run_provenance": "All three cases executed in one sequential invocation of this script",
        "reference": summarize(reference),
        "half_time_step": compare(reference, time_refined),
        "double_grid": compare(reference, mesh_refined),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if all(summary["budget_accepted"] for summary in (result["reference"], result["half_time_step"]["candidate"], result["double_grid"]["candidate"])) else 1


if __name__ == "__main__":
    raise SystemExit(main())
