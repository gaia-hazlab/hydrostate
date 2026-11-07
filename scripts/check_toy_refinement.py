"""Measure toy time/mesh sensitivity; this is not field validation."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from puget_modflow.toy_watershed import ToyConfig, run_toy


def compare(reference: dict, candidate: dict) -> dict:
    base = reference["records"]
    other = candidate["records"]
    if [r["time_h"] for r in base] != [r["time_h"] for r in other]:
        raise ValueError("Time windows must match")
    q0 = np.array([r["outlet_m3_s"] for r in base])
    q1 = np.array([r["outlet_m3_s"] for r in other])
    time = np.asarray([r["time_h"] for r in base])
    event = time >= 24
    theta0 = np.array([r["theta_root"] for r in base])
    theta1 = np.array([r["theta_root"] for r in other])
    # Compare pressure on the reference cell centers. The record's nearest-cell
    # pressure sensor changes support on refinement, so do not compare it directly.
    d0, x0 = np.meshgrid(reference["grid"]["depth_m"], reference["grid"]["x_m"], indexing="ij")
    query = np.column_stack([d0.ravel(), x0.ravel()])
    fine_psi = RegularGridInterpolator(
        (candidate["grid"]["depth_m"], candidate["grid"]["x_m"]),
        np.asarray(candidate["snapshots"][-1]["pressure_head_m"]),
        bounds_error=True,
    )(query).reshape(d0.shape)
    difference = fine_psi - np.asarray(reference["snapshots"][-1]["pressure_head_m"])
    return {
        "candidate_config": candidate["meta"]["config"],
        "budget_accepted": candidate["budget"]["accepted"],
        "budget_residual_mm": candidate["budget"]["signed_residual_mm"],
        "max_root_theta_difference_m3_m3": float(np.max(abs(theta1 - theta0))),
        "peak_outlet_relative_difference": float((q1.max() - q0.max()) / q0.max()),
        "peak_outlet_time_difference_h": float(other[int(q1.argmax())]["time_h"] - base[int(q0.argmax())]["time_h"]),
        "event_window_start_h": 24,
        "event_peak_outlet_relative_difference": float((q1[event].max() - q0[event].max()) / q0[event].max()) if event.any() and q0[event].max() > 0 else None,
        "event_peak_time_difference_h": float(time[event][int(q1[event].argmax())] - time[event][int(q0[event].argmax())]) if event.any() else None,
        "total_outlet_relative_difference": float((candidate["budget"]["outlet_m3"] - reference["budget"]["outlet_m3"]) / reference["budget"]["outlet_m3"]),
        "final_pressure_head_rmse_m": float(np.sqrt(np.mean(difference ** 2))),
        "final_pressure_head_max_difference_m": float(np.max(abs(difference))),
        "diagnostics": candidate["diagnostics"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("runs/toy_watershed/results.json"))
    parser.add_argument("--output", type=Path, default=Path("runs/toy_watershed/refinement.json"))
    args = parser.parse_args()
    baseline = json.loads(args.input.read_text())
    c = ToyConfig(**baseline["meta"]["config"])
    if baseline["meta"]["scenario"] != "baseline":
        raise ValueError("This benchmark compares the baseline scenario")
    result = {
        "interpretation": "Two sensitivity experiments, not formal mesh convergence or independent physical validation.",
        "reference_config": asdict(c),
        "reference_budget_accepted": baseline["budget"]["accepted"],
    }
    print("Running half-time-step comparison", flush=True)
    time_refined = run_toy(config=replace(c, max_step_s=c.max_step_s / 2), snapshot_interval_h=c.hours)
    result["half_time_step"] = compare(baseline, time_refined)
    print("Running doubled-horizontal-and-vertical-grid comparison", flush=True)
    space_refined = run_toy(config=replace(c, nx=c.nx * 2, nz=c.nz * 2), snapshot_interval_h=c.hours)
    result["double_grid"] = compare(baseline, space_refined)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
