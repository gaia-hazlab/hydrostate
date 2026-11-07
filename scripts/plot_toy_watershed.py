"""Export the synthetic Richards watershed run as figures and a CSV.

Run with the project environment; this does not invoke DHSVM or MODFLOW.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "runs" / ".matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def export_run(input_path: Path, output_dir: Path) -> dict[str, object]:
    data = json.loads(input_path.read_text(encoding="utf-8"))
    records = data["records"]
    snapshots = data["snapshots"]
    grid = data["grid"]
    output_dir.mkdir(parents=True, exist_ok=True)
    x = np.asarray(grid["x_m"])
    depth = np.asarray(grid["depth_m"])
    dx, dz = float(grid["dx_m"]), float(grid["dz_m"])
    xe = np.r_[x - dx / 2, x[-1] + dx / 2]
    de = np.r_[depth - dz / 2, depth[-1] + dz / 2]
    time = np.asarray([r["time_h"] for r in records])
    outflow = np.asarray([r["outlet_m3_s"] for r in records])
    peak = int(np.argmax(outflow))
    peak_t = float(time[peak])
    wettest_t = float(time[int(np.argmax([r["theta_root"] for r in records]))])
    selected = [snapshots[0], min(snapshots, key=lambda s: abs(s["time_h"] - wettest_t))]
    all_theta = np.asarray([s["theta"] for s in selected])
    low, high = float(all_theta.min()), float(all_theta.max())
    lens = np.asarray([["lens" in str(v).lower() for v in row] for row in grid["layer_name"]])

    fig = plt.figure(figsize=(13.2, 9), layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 0.95, 0.95])
    sections = [fig.add_subplot(gs[0, i]) for i in range(2)]
    for ax, snap in zip(sections, selected):
        mesh = ax.pcolormesh(xe, de, snap["theta"], cmap="cividis", vmin=low, vmax=high, shading="flat")
        pressure = np.asarray(snap["pressure_head_m"])
        if pressure.min() < 0 < pressure.max():
            ax.contour(x, depth, pressure, levels=[0], colors="#f2f2f2", linewidths=1.7)
            ax.contour(x, depth, pressure, levels=[0], colors="#1f2937", linewidths=0.7, linestyles="--")
        if lens.any() and not lens.all():
            ax.contourf(x, depth, lens.astype(float), levels=[0.5, 1.5], colors="none", hatches=["////"])
        ax.set(xlabel="Distance downslope (m)", ylabel="Depth normal to surface (m)",
               title=f"{'Initial state' if snap is selected[0] else 'Wettest root zone'} · hour {snap['time_h']:g}")
        ax.invert_yaxis()
    fig.colorbar(mesh, ax=sections, label="Volumetric water content (m³/m³)", shrink=0.9)
    sections[0].legend(handles=[Line2D([], [], color="#1f2937", ls="--", label="Zero pressure contour")],
                       loc="lower left", fontsize=9, framealpha=0.95)
    if lens.any():
        iz, ix = np.where(lens)
        sections[0].text(float(x[ix].mean()), float(depth[iz].mean()), "Low-K lens",
                         ha="center", va="center", fontsize=9,
                         bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85})

    ax = fig.add_subplot(gs[1, :])
    water_input = np.asarray([r["rain_mm_h"] + r["melt_mm_h"] for r in records])
    ax.fill_between(time, water_input, color="#8fb8ce", alpha=0.65, step="pre", label="Rain + melt")
    ax.set(ylabel="Rain + melt (mm/h)", xlim=(time[0], time[-1]))
    aq = ax.twinx()
    aq.plot(time, outflow * 1000, color="#1e5288", lw=1.8, drawstyle="steps-pre", label="Channel outlet (hourly mean)")
    aq.set_ylabel("Channel outlet (L/s)")
    ax.axvline(wettest_t, color="#555555", ls=":", lw=1)
    ax.legend(loc="upper left", frameon=False)
    aq.legend(loc="upper right", frameon=False)
    ax.grid(axis="y", alpha=0.2)

    ax = fig.add_subplot(gs[2, 0])
    ax.plot(time, [r["theta_root"] for r in records], color="#39794b", lw=1.8)
    ax.set(xlabel="Elapsed time (h)", ylabel="Root-zone water content\n(m³/m³)", xlim=(time[0], time[-1]))
    ax.grid(alpha=0.2)
    ax = fig.add_subplot(gs[2, 1])
    ax.plot(time, [r["pressure_sensor_kpa"] for r in records], color="#aa5a26", lw=1.8)
    ax.axhline(0, color="#555555", lw=0.8)
    ax.set(xlabel="Elapsed time (h)", ylabel="Shallow pore pressure (kPa)", xlim=(time[0], time[-1]))
    ax.grid(alpha=0.2)
    fig.suptitle("Synthetic headwater hillslope and channel · 2-D mixed Richards prototype", fontsize=15)
    fig.savefig(output_dir / "summary.png", dpi=170, bbox_inches="tight", pad_inches=0.15)
    fig.savefig(output_dir / "summary.pdf", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    with (output_dir / "timeseries.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    summary = {
        "model": data["meta"],
        "budget": data["budget"],
        "diagnostics": data["diagnostics"],
        "peak_outlet_l_s": float(outflow.max() * 1000),
        "peak_outlet_time_h": peak_t,
        "wettest_root_zone_time_h": wettest_t,
        "event_peak_outlet_l_s": float(outflow[time >= 24].max() * 1000) if (time >= 24).any() else None,
        "event_peak_outlet_time_h": float(time[time >= 24][int(np.argmax(outflow[time >= 24]))]) if (time >= 24).any() else None,
        "event_window_start_h": 24,
        "root_theta_range": [float(min(r["theta_root"] for r in records)),
                             float(max(r["theta_root"] for r in records))],
        "sensor_pore_pressure_range_kpa": [float(min(r["pressure_sensor_kpa"] for r in records)),
                                           float(max(r["pressure_sensor_kpa"] for r in records))],
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("runs/toy_watershed/results.json"))
    parser.add_argument("--output", type=Path, default=Path("runs/toy_watershed"))
    args = parser.parse_args()
    summary = export_run(args.input, args.output)
    print(json.dumps({k: v for k, v in summary.items() if k not in {"model", "budget", "diagnostics"}}, indent=2))
    print(f"Figures and CSV: {args.output.resolve()}")


if __name__ == "__main__":
    main()
