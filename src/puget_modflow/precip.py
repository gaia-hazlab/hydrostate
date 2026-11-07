from __future__ import annotations

from datetime import date
from pathlib import Path
import csv
import re
import zipfile


MM_PER_INCH = 25.4


def precip_mm_to_inches(value: float) -> float:
    return value / MM_PER_INCH


def prism_date_from_name(name: str) -> date:
    """Extract a monthly or daily PRISM date token from a file name."""
    matches = re.findall(r"(?:^|_)(\d{6}|\d{8})(?:_|\.|$)", Path(name).name)
    if not matches:
        raise ValueError(f"Could not find PRISM date token in {name}")
    token = matches[-1]
    if len(token) == 6:
        return date(int(token[:4]), int(token[4:6]), 1)
    return date(int(token[:4]), int(token[4:6]), int(token[6:8]))


def find_prism_grids(input_dir: str | Path) -> list[Path]:
    root = Path(input_dir)
    if not root.exists():
        raise FileNotFoundError(root)
    candidates = sorted(root.glob("*.zip")) + sorted(root.glob("*.tif"))
    return [path for path in candidates if path.is_file()]


def aggregate_prism_to_domain_csv(
    *,
    input_dir: str | Path,
    output_csv: str | Path,
    bounds_wgs84: tuple[float, float, float, float],
) -> Path:
    """Aggregate PRISM COG or zipped COG grids to domain-mean precipitation.

    Output units include both PRISM-native millimeters and inches for the
    recharge conversion workflow.
    """
    rasterio, from_bounds, transform_bounds = _require_rasterio()

    grids = find_prism_grids(input_dir)
    if not grids:
        raise FileNotFoundError(f"No .zip or .tif PRISM grids found in {input_dir}")

    rows: list[dict[str, str]] = []
    for grid in grids:
        raster_path = _rasterio_path(grid)
        with rasterio.open(raster_path) as dataset:
            if dataset.crs is None:
                raise ValueError(f"Raster has no CRS: {grid}")
            left, bottom, right, top = transform_bounds(
                "EPSG:4326",
                dataset.crs,
                bounds_wgs84[0],
                bounds_wgs84[1],
                bounds_wgs84[2],
                bounds_wgs84[3],
                densify_pts=21,
            )
            window = from_bounds(left, bottom, right, top, dataset.transform)
            data = dataset.read(1, window=window, masked=True)
            if data.count() == 0:
                raise ValueError(f"No raster cells overlap configured domain for {grid}")
            precip_mm = float(data.mean())
        dt = prism_date_from_name(grid.name)
        rows.append(
            {
                "date": dt.isoformat(),
                "source_file": str(grid),
                "precip_mm": f"{precip_mm:.12g}",
                "precip_in": f"{precip_mm_to_inches(precip_mm):.12g}",
            }
        )

    output = Path(output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["date", "source_file", "precip_mm", "precip_in"],
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["date"]))
    return output


def _rasterio_path(path: Path) -> str:
    if path.suffix.lower() == ".tif":
        return str(path)
    if path.suffix.lower() != ".zip":
        raise ValueError(f"Unsupported PRISM raster archive: {path}")
    with zipfile.ZipFile(path) as archive:
        tif_names = [name for name in archive.namelist() if name.lower().endswith(".tif")]
    if not tif_names:
        raise ValueError(f"No .tif found inside {path}")
    return f"zip://{path.resolve()}!{tif_names[0]}"


def _require_rasterio():
    try:
        import rasterio
        from rasterio.windows import from_bounds
        from rasterio.warp import transform_bounds
    except ImportError as exc:
        raise RuntimeError(
            "Raster aggregation requires rasterio. Install with "
            '`python -m pip install -e ".[geo]"`.'
        ) from exc
    return rasterio, from_bounds, transform_bounds
