"""Render the synthetic valley and its observation comparisons with Matplotlib."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("runs/valley3d/results.json"))
    parser.add_argument("--sensors", type=Path, default=Path("examples/sensors/sensors.csv"))
    parser.add_argument("--comparison", type=Path, default=Path("runs/sensors/report.json"))
    parser.add_argument("--output", type=Path, default=Path("runs/valley3d/valley3d.png"))
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    grid = data["grid"]
    x, y, z = (np.array(grid[f"{axis}_m"]) for axis in "xyz")
    dx, dy, dz = (grid[f"d{axis}_m"] for axis in "xyz")
    surface = np.array(grid["surface_elevation_m"])
    final = data["snapshots"][-1]
    psi = np.array(final["pressure_head_m"], dtype=float)
    theta = np.array(final["theta"], dtype=float)
    layers = np.array(grid["layer_name"])
    records = data["records"]
    hours = np.array([r["time_h"] for r in records])
    sensors = []
    if args.sensors.exists():
        with args.sensors.open(newline="") as stream:
            sensors = list(csv.DictReader(stream))
    comparison = json.loads(args.comparison.read_text()) if args.comparison.exists() else None
    plt.rcParams.update({"font.size": 10, "axes.titleweight": "bold", "axes.spines.top": False})
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), layout="constrained")
    fig.suptitle("Hydrostate · synthetic 3D valley", fontsize=20, fontweight="bold")

    ax = axes[0, 0]
    image = ax.pcolormesh(x, y, surface, shading="nearest", cmap="terrain", vmin=3.0, vmax=6.0)
    fig.colorbar(image, ax=ax, label="Surface elevation (m)", shrink=.8)
    locations = {}
    for row in sensors:
        point = float(row["x_m"]), float(row["y_m"])
        locations.setdefault(point, []).append(row)
    for (sx, sy), rows in locations.items():
        ax.scatter(sx, sy, marker="o" if rows[0]["role"] == "validation" else "s", s=35,
                   facecolor="white", edgecolor="#122e3e", zorder=3)
        ax.annotate(" / ".join(row["sensor_id"] for row in rows), (sx, sy),
                    xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set(xlabel="Eastward x (m)", ylabel="Cross-valley y (m)", title="A  Stepped surface and synthetic sensors",
           xlim=(x[0] - dx / 2, x[-1] + dx / 2 + 10), ylim=(y[0] - dy / 2, y[-1] + dy / 2))

    j, i = len(y) // 2, len(x) // 2
    ax = axes[0, 1]
    image = ax.pcolormesh(x, z, psi[:, j, :], shading="nearest", cmap="RdBu", vmin=-3.5, vmax=3.5)
    ax.contour(x, z, psi[:, j, :], levels=[0], colors="black", linewidths=1.2)
    lens = (layers[:, j, :] == "low_K_lens").astype(float)
    if lens.any():
        ax.contour(x, z, lens, levels=[.5], colors="#c18400", linewidths=2)
    ax.step(np.r_[x - dx / 2, x[-1] + dx / 2], np.r_[surface[j], surface[j, -1]], where="post", color="#333333", lw=1)
    fig.colorbar(image, ax=ax, label="Pressure head ψ (m)", shrink=.8)
    ax.set(xlabel="Eastward x (m)", ylabel="Elevation z (m)", title=f"B  Pressure at {hours[-1]:g} h; y = {y[j]:g} m",
           ylim=(z[0] - dz / 2, z[-1] + dz / 2))
    ax.text(.03, .95, "Black: ψ = 0\nGold: low-K lens", transform=ax.transAxes, va="top", fontsize=8)

    ax = axes[0, 2]
    image = ax.pcolormesh(y, z, theta[:, :, i], shading="nearest", cmap="YlGnBu", vmin=.2, vmax=.46)
    ax.contour(y, z, psi[:, :, i], levels=[0], colors="black", linewidths=1.2)
    fig.colorbar(image, ax=ax, label="Water content (m³/m³)", shrink=.8)
    ax.set(xlabel="Cross-valley y (m)", ylabel="Elevation z (m)", title=f"C  Moisture at {hours[-1]:g} h; x = {x[i]:g} m",
           ylim=(z[0] - dz / 2, z[-1] + dz / 2))

    ax = axes[1, 0]
    ax.bar(hours[1:] - .5, [r["precip_mm_h"] for r in records[1:]], width=1,
           color="#a2cbe0", label="Infiltrating liquid rain")
    ax.set(xlabel="Elapsed time (h)", ylabel="Liquid rain (mm/h)", title="D  Forcing and open-boundary response", ylim=(0, .25))
    right = ax.twinx()
    right.plot(hours[1:], [1000 * r["outlet_m3_s"] for r in records[1:]], color="#12616a", lw=2,
               label="Net east-boundary flux")
    right.set_ylabel("Net east-boundary outflow (L/s)", color="#12616a")
    right.spines["top"].set_visible(False)
    ax.text(.03, .95, "Boundary flux is a groundwater surrogate", transform=ax.transAxes, va="top", fontsize=8)

    ax = axes[1, 1]
    if comparison:
        for identifier, color in (("SM_1", "#926b12"), ("SM_2", "#356c96")):
            matches = [r for r in comparison["matches"] if r["sensor_id"] == identifier]
            if not matches:
                continue
            from datetime import datetime
            start = datetime.fromisoformat(data["meta"]["start_utc"].replace("Z", "+00:00"))
            t = [(datetime.fromisoformat(r["valid_time_utc"].replace("Z", "+00:00")) - start).total_seconds() / 3600 for r in matches]
            ax.plot(t, [r["predicted"] for r in matches], color=color, label=f"{identifier} model")
            ax.errorbar(t, [r["observed"] for r in matches], yerr=[r["uncertainty"] for r in matches],
                        fmt="o", ms=4, color=color, capsize=3, label=f"{identifier} synthetic")
        ax.legend(fontsize=8, loc="best")
    else:
        ax.plot(hours, [r["theta_surface"] for r in records], color="#356c96", label="Mean top-cell moisture")
        ax.legend(fontsize=8)
    ax.set(xlabel="Elapsed time (h)", ylabel="Water content (m³/m³)", title="E  Observation operators and held-out site")
    ax.text(.03, .03, "SM_2 is held out; offsets are deterministic", transform=ax.transAxes, fontsize=8)

    ax = axes[1, 2]
    ax.plot(hours, [r["budget_residual_mm"] for r in records], color="#12616a", lw=2)
    ax.axhline(0, color="gray", lw=.6)
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.set(xlabel="Elapsed time (h)", ylabel="Cumulative balance residual (mm)", title="F  Conservation audit")
    d = data["diagnostics"]
    ax.text(.04, .04, f"{d['active_cells']} active Cartesian cells\n"
            f"Max local residual: {d['max_local_residual_m3_s']:.2g} m³/s\n"
            "Valley mesh convergence remains unestablished", transform=ax.transAxes, va="bottom", fontsize=9)
    fig.supxlabel("Hypothetical local coordinates and parameters · no field validation, assimilation or hazard prediction", fontsize=10)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
