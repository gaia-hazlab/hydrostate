#!/usr/bin/env python3
"""Embed compact, precomputed toy outputs into the existing inline viewer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re


def compact_data(source: dict, stride: int) -> dict:
    """Keep the full grid, hourly traces, and a rounded subset of frames."""
    grid = source["grid"]
    names = grid["layer_name"]
    result_grid = {
        key: grid[key]
        for key in (
            "x_m", "depth_m", "surface_elevation_m", "dx_m", "dz_m",
            "width_m", "slope_rise_run",
        )
    }
    result_grid["lens_mask"] = [
        ["lens" in str(name).lower() for name in row] for row in names
    ]
    snapshots = source["snapshots"]
    chosen = list(range(0, len(snapshots), stride))
    if chosen[-1] != len(snapshots) - 1:
        chosen.append(len(snapshots) - 1)
    frames = [
        {
            "time_h": round(snapshots[i]["time_h"], 4),
            "theta": [[round(v, 4) for v in row] for row in snapshots[i]["theta"]],
            "pressure_head_m": [
                [round(v, 3) for v in row]
                for row in snapshots[i]["pressure_head_m"]
            ],
        }
        for i in chosen
    ]
    keep = (
        "time_h", "rain_mm_h", "melt_mm_h", "theta_root", "outlet_m3_s",
        "root_stress",
    )
    records = [
        {key: round(row[key], 8) for key in keep} for row in source["records"]
    ]
    wettest = max(records, key=lambda row: row["theta_root"])["time_h"]
    return {
        "grid": result_grid,
        "records": records,
        "snapshots": frames,
        "default_frame": min(
            range(len(frames)), key=lambda i: abs(frames[i]["time_h"] - wettest)
        ),
        "meta": {
            "label": "Precomputed synthetic Richards watershed toy; not DHSVM/MODFLOW",
            "frame_stride": stride,
            "field_rounding": {"theta_decimals": 4, "pressure_head_m_decimals": 3},
            "time_definition": source.get("meta", {}).get("time_definition"),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--fragment", type=Path, required=True,
                        help="Existing visualization fragment whose embedded data will be updated")
    parser.add_argument("--stride", type=int, default=2)
    parser.add_argument("--max-bytes", type=int, default=950_000)
    args = parser.parse_args()
    source = json.loads(args.results.read_text())
    if not source.get("snapshots") or not source.get("records"):
        raise ValueError("Numerical output must include snapshots and records")
    fragment = args.fragment.read_text()
    pattern = re.compile(
        r'(?<=<script type="application/json" id="toy-numerical-data">).*?(?=</script>)',
        flags=re.DOTALL,
    )
    if len(pattern.findall(fragment)) != 1:
        raise ValueError("Expected one inline numerical-data placeholder")
    stride = max(1, args.stride)
    while True:
        data = compact_data(source, stride)
        encoded = json.dumps(data, separators=(",", ":"), allow_nan=False)
        encoded = encoded.replace("<", "\\u003c")
        output = pattern.sub(lambda _: encoded, fragment)
        size = len(output.encode())
        if size <= args.max_bytes:
            args.fragment.write_text(output)
            print(json.dumps({"fragment": str(args.fragment), "bytes": size,
                              "frames": len(data["snapshots"]), "stride": stride,
                              "hourly_records": len(data["records"])}))
            return
        if stride >= len(source["snapshots"]):
            raise ValueError("Even endpoint frames exceed the viewer size budget")
        stride += 1


if __name__ == "__main__":
    main()
