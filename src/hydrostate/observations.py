"""CSV sensor ingestion and support-aware forward comparisons; no assimilation.

Coordinates are Cartesian metres with z positive upward. This module uses the
standard library so that CSV QA and analytical operator tests stay lightweight.
"""

from __future__ import annotations

from bisect import bisect_left
import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import json
import hashlib
import math
from pathlib import Path
import sys
from typing import Any, Iterable


KINDS = {
    "soil_moisture", "pressure_head", "pore_pressure", "hydraulic_head",
    "water_level_depth", "streamflow", "precipitation",
}
RATE_KINDS = {"streamflow", "precipitation"}
HEAD_KINDS = {"hydraulic_head", "water_level_depth"}
CANONICAL_UNITS = {
    "soil_moisture": "m3/m3", "pressure_head": "m", "pore_pressure": "Pa",
    "hydraulic_head": "m", "water_level_depth": "m", "streamflow": "m3/s",
    "precipitation": "mm/h",
}
QUALITY_FLAGS = {"good", "provisional", "suspect", "bad", "missing"}
SENSOR_COLUMNS = [
    "sensor_id", "kind", "x_m", "y_m", "z_m", "depth_m", "screen_top_m",
    "screen_bottom_m", "support_x_m", "support_y_m", "support_z_m", "units",
    "uncertainty", "horizontal_crs", "vertical_datum", "aquifer_id",
    "screen_assumption", "outlet_id", "role", "provenance",
]
OBSERVATION_COLUMNS = [
    "sensor_id", "valid_time_utc", "arrival_time_utc", "interval_start_utc",
    "value", "units", "uncertainty", "quality_flag",
]


class ObservationError(ValueError):
    """Malformed observations or an undefined/mismatched forward operator."""


def utc_time(value: str | datetime) -> datetime:
    """Parse an explicit timezone; never guess the timezone of a naive date."""
    try:
        result = value if isinstance(value, datetime) else datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (TypeError, ValueError) as exc:
        raise ObservationError(f"Invalid timestamp: {value!r}") from exc
    if result.tzinfo is None or result.utcoffset() is None:
        raise ObservationError(f"Timestamp needs an explicit timezone: {value!r}")
    return result.astimezone(timezone.utc)


def iso_utc(value: datetime) -> str:
    return utc_time(value).isoformat().replace("+00:00", "Z")


def _finite(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ObservationError(f"{label} must be numeric") from exc
    if not math.isfinite(result):
        raise ObservationError(f"{label} must be finite")
    return result


def _optional(row: dict[str, str], key: str) -> float | None:
    return _finite(row[key], key) if row.get(key, "").strip() else None


def unit_factor(kind: str, units: str, interval_seconds: float | None = None) -> float:
    """Multiplicative conversion into the canonical units; no datum conversion."""
    aliases = {
        "m³/m³": "m3/m3", "m^3/m^3": "m3/m3", "fraction": "m3/m3",
        "m³/s": "m3/s", "m^3/s": "m3/s", "l/s": "L/s",
        "ft3/s": "cfs", "ft^3/s": "cfs", "ft³/s": "cfs",
        "percent": "%", "kpa": "kPa", "pa": "Pa",
    }
    units = aliases.get(units.strip(), units.strip())
    if kind == "soil_moisture":
        choices = {"m3/m3": 1.0, "%": 0.01}
    elif kind in {"pressure_head", "hydraulic_head", "water_level_depth"}:
        choices = {"m": 1.0, "cm": 0.01, "mm": 0.001, "ft": 0.3048}
    elif kind == "pore_pressure":
        choices = {"Pa": 1.0, "kPa": 1000.0}
    elif kind == "streamflow":
        choices = {"m3/s": 1.0, "L/s": 0.001, "cfs": 0.028316846592}
    elif kind == "precipitation":
        choices = {"mm/h": 1.0, "m/s": 3_600_000.0, "in/h": 25.4}
        if units in {"mm", "in"}:
            if interval_seconds is None or interval_seconds <= 0:
                raise ObservationError("Precipitation accumulation needs a positive observation interval")
            return (1.0 if units == "mm" else 25.4) * 3600.0 / interval_seconds
    else:
        raise ObservationError(f"Unsupported sensor kind: {kind}")
    if units not in choices:
        raise ObservationError(f"Unsupported units {units!r} for {kind}")
    return choices[units]


@dataclass(frozen=True)
class Sensor:
    sensor_id: str
    kind: str
    x_m: float
    y_m: float
    units: str
    uncertainty: float
    horizontal_crs: str
    vertical_datum: str
    role: str
    z_m: float | None = None
    depth_m: float | None = None
    screen_top_m: float | None = None
    screen_bottom_m: float | None = None
    support_x_m: float = 0.0
    support_y_m: float = 0.0
    support_z_m: float = 0.0
    aquifer_id: str = ""
    screen_assumption: str = ""
    outlet_id: str = ""
    provenance: str = ""


@dataclass(frozen=True)
class Observation:
    sensor_id: str
    valid_time: datetime
    arrival_time: datetime
    value: float | None
    uncertainty: float
    units: str
    quality_flag: str
    interval_start: datetime | None = None


def _rows(path: str | Path, required: Iterable[str]) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames and len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ObservationError(f"{path}: duplicate CSV header names")
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ObservationError(f"{path}: missing columns {sorted(set(required) - set(reader.fieldnames or []))}")
        result = []
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ObservationError(f"{path}:{line}: inconsistent CSV column count")
            result.append({key: value.strip() for key, value in row.items()})
        return result


def load_sensors(path: str | Path) -> dict[str, Sensor]:
    required = {"sensor_id", "kind", "x_m", "y_m", "units", "uncertainty", "horizontal_crs", "vertical_datum", "role"}
    sensors: dict[str, Sensor] = {}
    for row in _rows(path, required):
        identifier = row["sensor_id"]
        if not identifier or identifier in sensors:
            raise ObservationError(f"Empty or duplicate sensor_id: {identifier!r}")
        kind = row["kind"]
        if kind not in KINDS:
            raise ObservationError(f"{identifier}: unsupported kind {kind!r}")
        if row["role"] not in {"train", "validation"}:
            raise ObservationError(f"{identifier}: role must be train or validation")
        if not row["horizontal_crs"] or not row["vertical_datum"]:
            raise ObservationError(f"{identifier}: CRS and datum are required")
        sensor = Sensor(
            sensor_id=identifier, kind=kind, x_m=_finite(row["x_m"], "x_m"),
            y_m=_finite(row["y_m"], "y_m"), units=row["units"],
            uncertainty=_finite(row["uncertainty"], "uncertainty"),
            horizontal_crs=row["horizontal_crs"], vertical_datum=row["vertical_datum"], role=row["role"],
            **{key: _optional(row, key) for key in ("z_m", "depth_m", "screen_top_m", "screen_bottom_m")},
            **{key: _optional(row, key) or 0.0 for key in ("support_x_m", "support_y_m", "support_z_m")},
            **{key: row.get(key, "") for key in ("aquifer_id", "screen_assumption", "outlet_id", "provenance")},
        )
        if sensor.uncertainty <= 0:
            raise ObservationError(f"{identifier}: uncertainty must be positive")
        if min(sensor.support_x_m, sensor.support_y_m, sensor.support_z_m) < 0:
            raise ObservationError(f"{identifier}: support widths cannot be negative")
        has_screen = sensor.screen_top_m is not None or sensor.screen_bottom_m is not None
        if has_screen:
            if kind not in HEAD_KINDS or sensor.screen_top_m is None or sensor.screen_bottom_m is None:
                raise ObservationError(f"{identifier}: a complete screen interval is only supported for head/level wells")
            if sensor.screen_top_m <= sensor.screen_bottom_m:
                raise ObservationError(f"{identifier}: screen top must be above screen bottom")
            if not sensor.aquifer_id or sensor.screen_assumption != "single_aquifer_uniform_conductance":
                raise ObservationError(f"{identifier}: declare the screen's aquifer and single_aquifer_uniform_conductance assumption")
            if sensor.z_m is not None or sensor.depth_m is not None or sensor.support_z_m != 0:
                raise ObservationError(f"{identifier}: use the screen interval instead of z/depth/vertical support")
        elif kind not in RATE_KINDS:
            if (sensor.z_m is None) == (sensor.depth_m is None):
                raise ObservationError(f"{identifier}: provide exactly one z_m or depth_m")
        if sensor.depth_m is not None and sensor.depth_m < 0:
            raise ObservationError(f"{identifier}: depth_m must be positive below the local surface")
        if kind == "streamflow" and not sensor.outlet_id:
            raise ObservationError(f"{identifier}: streamflow needs a model outlet_id")
        unit_factor(kind, sensor.units, 3600.0)
        sensors[identifier] = sensor
    if not sensors:
        raise ObservationError("Sensor metadata is empty")
    return sensors


def load_observations(path: str | Path, sensors: dict[str, Sensor]) -> list[Observation]:
    required = {"sensor_id", "valid_time_utc", "arrival_time_utc", "value", "quality_flag"}
    observations = []
    seen: set[tuple[str, datetime]] = set()
    for row in _rows(path, required):
        if row["sensor_id"] not in sensors:
            raise ObservationError(f"Unknown sensor_id: {row['sensor_id']!r}")
        sensor = sensors[row["sensor_id"]]
        valid, arrival = utc_time(row["valid_time_utc"]), utc_time(row["arrival_time_utc"])
        if arrival < valid:
            raise ObservationError(f"{sensor.sensor_id}: arrival precedes valid time")
        duplicate_key = (sensor.sensor_id, valid)
        if duplicate_key in seen:
            raise ObservationError(f"Duplicate sensor/time: {sensor.sensor_id}, {iso_utc(valid)}")
        seen.add(duplicate_key)
        start = utc_time(row["interval_start_utc"]) if row.get("interval_start_utc") else None
        if sensor.kind in RATE_KINDS:
            if start is None or start >= valid:
                raise ObservationError(f"{sensor.sensor_id}: rates/accumulations need interval_start_utc < valid_time_utc")
        elif start is not None:
            raise ObservationError(f"{sensor.sensor_id}: temporal means of state observations are not implemented")
        duration = (valid - start).total_seconds() if start is not None else None
        units = row.get("units") or sensor.units
        factor = unit_factor(sensor.kind, units, duration)
        quality = row["quality_flag"]
        if quality not in QUALITY_FLAGS:
            raise ObservationError(f"{sensor.sensor_id}: unrecognized quality flag {quality!r}")
        raw_value = row["value"]
        missing = raw_value.lower() in {"", "na", "nan", "null", "none"}
        value = None if missing else _finite(raw_value, "observation value") * factor
        if quality == "missing" and value is not None:
            raise ObservationError(f"{sensor.sensor_id}: missing quality flag conflicts with a numeric value")
        if value is not None and sensor.kind == "soil_moisture" and not 0 <= value <= 1:
            raise ObservationError(f"{sensor.sensor_id}: volumetric water content must lie in [0, 1]")
        if value is not None and sensor.kind == "precipitation" and value < 0:
            raise ObservationError(f"{sensor.sensor_id}: precipitation cannot be negative")
        sigma = (_finite(row["uncertainty"], "uncertainty") * factor if row.get("uncertainty")
                 else sensor.uncertainty * unit_factor(sensor.kind, sensor.units, duration))
        if sigma <= 0:
            raise ObservationError(f"{sensor.sensor_id}: uncertainty must be positive")
        observations.append(Observation(sensor.sensor_id, valid, arrival, value, sigma,
                                        CANONICAL_UNITS[sensor.kind], quality, start))
    return sorted(observations, key=lambda o: (o.valid_time, o.sensor_id))


def _axis_weights(axis: list[float], spacing: float, center: float, width: float) -> list[tuple[int, float]]:
    low, high = axis[0] - spacing / 2, axis[-1] + spacing / 2
    tolerance = max(1.0, abs(low), abs(high)) * 1e-10
    if center - width / 2 < low - tolerance or center + width / 2 > high + tolerance:
        raise ObservationError("Sensor support lies outside the model grid")
    if width > 0:
        overlaps = [(i, max(0.0, min(center + width / 2, x + spacing / 2)
                                 - max(center - width / 2, x - spacing / 2))) for i, x in enumerate(axis)]
        return [(i, overlap / width) for i, overlap in overlaps if overlap > 1e-12]
    if center <= axis[0]:
        return [(0, 1.0)]
    if center >= axis[-1]:
        return [(len(axis) - 1, 1.0)]
    upper = bisect_left(axis, center)
    fraction = (center - axis[upper - 1]) / (axis[upper] - axis[upper - 1])
    return [(i, weight) for i, weight in [(upper - 1, 1 - fraction), (upper, fraction)] if weight > 1e-12]


class ModelSampler:
    """Forward operators for the documented Cartesian 3-D JSON result schema."""

    def __init__(self, model: dict[str, Any]):
        self.model, self.meta, self.grid = model, model["meta"], model["grid"]
        if "start_utc" not in self.meta:
            raise ObservationError("Model meta.start_utc is required")
        self.start = utc_time(self.meta["start_utc"])
        self.axes = [list(map(float, self.grid[f"{name}_m"])) for name in ("x", "y", "z")]
        self.spacing = [_finite(self.grid[f"d{name}_m"], f"d{name}_m") for name in ("x", "y", "z")]
        for axis, spacing in zip(self.axes, self.spacing):
            if not axis or spacing <= 0 or any(not math.isfinite(v) for v in axis):
                raise ObservationError("Grid axes must be finite and spacings positive")
            if any(b <= a or not math.isclose(b - a, spacing, rel_tol=1e-8, abs_tol=1e-9) for a, b in zip(axis, axis[1:])):
                raise ObservationError("Only ascending, uniformly spaced Cartesian axes are supported")
        self.snapshots = sorted(model["snapshots"], key=lambda row: row["time_h"])
        self.records = sorted(model["records"], key=lambda row: row["time_h"])
        self.times = [float(row["time_h"]) for row in self.snapshots]
        self.record_times = [float(row["time_h"]) for row in self.records]
        for times in (self.times, self.record_times):
            if not times or any(not math.isfinite(t) for t in times) or any(b <= a for a, b in zip(times, times[1:])):
                raise ObservationError("Model times must be unique, finite and increasing")
        nx, ny, nz = map(len, self.axes)
        active = self.grid["active"]
        if len(active) != nz or any(len(plane) != ny or any(len(row) != nx for row in plane) for plane in active):
            raise ObservationError("Active mask must have shape (nz, ny, nx)")

    @classmethod
    def from_file(cls, path: str | Path) -> ModelSampler:
        return cls(json.loads(Path(path).read_text()))

    def _reference(self, sensor: Sensor) -> None:
        reference = self.meta.get("coordinate_reference", {})
        if reference.get("horizontal_crs") != sensor.horizontal_crs:
            raise ObservationError(f"{sensor.sensor_id}: horizontal CRS mismatch; transform explicitly before ingestion")
        if reference.get("vertical_datum") != sensor.vertical_datum:
            raise ObservationError(f"{sensor.sensor_id}: vertical datum mismatch; transform explicitly before ingestion")
        if sensor.screen_top_m is not None and self.meta.get("aquifer_id") != sensor.aquifer_id:
            raise ObservationError(f"{sensor.sensor_id}: screen aquifer does not match the model's declared single aquifer")

    def _local_surface(self, sensor: Sensor) -> float:
        indices = []
        for axis, spacing, location in zip(self.axes[:2], self.spacing[:2], (sensor.x_m, sensor.y_m)):
            _axis_weights(axis, spacing, location, 0.0)
            indices.append(min(len(axis) - 1, max(0, int(math.floor((location - axis[0] + spacing / 2) / spacing)))))
        i, j = indices
        return _finite(self.grid["surface_elevation_m"][j][i], "local surface elevation")

    def _support(self, sensor: Sensor) -> list[tuple[int, int, int, float]]:
        wx = _axis_weights(self.axes[0], self.spacing[0], sensor.x_m, sensor.support_x_m)
        wy = _axis_weights(self.axes[1], self.spacing[1], sensor.y_m, sensor.support_y_m)
        if sensor.screen_top_m is not None:
            assert sensor.screen_bottom_m is not None
            center = (sensor.screen_top_m + sensor.screen_bottom_m) / 2
            width = sensor.screen_top_m - sensor.screen_bottom_m
        else:
            center = sensor.z_m if sensor.z_m is not None else self._local_surface(sensor) - float(sensor.depth_m)
            width = sensor.support_z_m
        wz = _axis_weights(self.axes[2], self.spacing[2], center, width)
        support = [(k, j, i, a * b * c) for i, a in wx for j, b in wy for k, c in wz]
        if any(not self.grid["active"][k][j][i] for k, j, i, _ in support):
            raise ObservationError(f"{sensor.sensor_id}: support intersects inactive/air cells; it is not renormalized")
        if not math.isclose(sum(w for *_, w in support), 1.0, abs_tol=1e-9):
            raise ObservationError(f"{sensor.sensor_id}: incomplete spatial support")
        return support

    def _state(self, sensor: Sensor, snapshot: dict[str, Any]) -> float:
        support = self._support(sensor)
        key = {"soil_moisture": "theta", "pressure_head": "pressure_head_m",
               "pore_pressure": "pressure_head_m", "hydraulic_head": "total_head_m",
               "water_level_depth": "total_head_m"}[sensor.kind]
        if sensor.screen_top_m is not None:
            for k, j, i, _ in support:
                if _finite(snapshot["pressure_head_m"][k][j][i], "screen pressure") < 0:
                    raise ObservationError(f"{sensor.sensor_id}: partially unsaturated screen is unsupported")
        try:
            result = sum(weight * _finite(snapshot[key][k][j][i], key) for k, j, i, weight in support)
        except (KeyError, IndexError, TypeError) as exc:
            raise ObservationError(f"Missing/malformed model field {key}") from exc
        if sensor.kind == "water_level_depth":
            return self._local_surface(sensor) - result
        if sensor.kind == "pore_pressure":
            return result * 1000.0 * 9.80665
        if sensor.kind == "soil_moisture" and not 0 <= result <= 1:
            raise ObservationError("Modeled volumetric water content lies outside [0, 1]")
        return result

    def predict(self, sensor: Sensor, observation: Observation) -> float:
        self._reference(sensor)
        end = (observation.valid_time - self.start).total_seconds() / 3600
        if sensor.kind in RATE_KINDS:
            if observation.interval_start is None:
                raise ObservationError("Rate observation requires an explicit averaging interval")
            start = (observation.interval_start - self.start).total_seconds() / 3600
            for axis, spacing, location in zip(self.axes[:2], self.spacing[:2], (sensor.x_m, sensor.y_m)):
                _axis_weights(axis, spacing, location, 0.0)
            if start < self.record_times[0] or end > self.record_times[-1] or start >= end:
                raise ObservationError("Observation interval lies outside modeled rate coverage")
            if sensor.kind == "streamflow":
                if not sensor.outlet_id or sensor.outlet_id != self.meta.get("outlet", {}).get("id"):
                    raise ObservationError(f"{sensor.sensor_id}: outlet_id does not match the modeled discharge")
                key = "outlet_m3_s"
            else:
                if self.meta.get("precipitation_spatial_support") != "uniform":
                    raise ObservationError("A domain-mean precipitation record is not a point-gauge operator without uniform forcing")
                key = "precip_mm_h"
            integral = coverage = 0.0
            for left, right in zip(self.records, self.records[1:]):
                overlap = max(0.0, min(end, right["time_h"]) - max(start, left["time_h"]))
                if overlap:
                    integral += overlap * _finite(right[key], key)
                    coverage += overlap
            if not math.isclose(coverage, end - start, abs_tol=1e-9):
                raise ObservationError("Incomplete rate interval coverage")
            return integral / (end - start)
        if end < self.times[0] or end > self.times[-1]:
            raise ObservationError("Observation time lies outside model snapshots")
        upper = bisect_left(self.times, end)
        if self.times[upper] == end:
            return self._state(sensor, self.snapshots[upper])
        lower = upper - 1
        fraction = (end - self.times[lower]) / (self.times[upper] - self.times[lower])
        return (1 - fraction) * self._state(sensor, self.snapshots[lower]) + fraction * self._state(sensor, self.snapshots[upper])


def metrics(rows: list[dict[str, Any]]) -> dict[str, float | int]:
    if not rows:
        return {"n": 0}
    residuals = [row["predicted"] - row["observed"] for row in rows]
    weights = [1 / row["uncertainty"] ** 2 for row in rows]
    return {
        "n": len(rows), "bias_model_minus_observed": sum(residuals) / len(rows),
        "rmse": math.sqrt(sum(r * r for r in residuals) / len(rows)),
        "mae": sum(abs(r) for r in residuals) / len(rows),
        "uncertainty_weighted_rmse": math.sqrt(sum(w * r * r for w, r in zip(weights, residuals)) / sum(weights)),
        "mean_squared_observation_normalized_residual": sum((r / row["uncertainty"]) ** 2 for r, row in zip(residuals, rows)) / len(rows),
        "fraction_within_observation_1_96sigma": sum(abs(r) <= 1.96 * row["uncertainty"] for r, row in zip(residuals, rows)) / len(rows),
    }


def compare_observations(
    model: dict[str, Any], sensors: dict[str, Sensor], observations: list[Observation],
    *, as_of: datetime | str | None = None, include_quality: tuple[str, ...] = ("good",),
) -> dict[str, Any]:
    """Compare forward predictions; holdouts and precipitation remain separate."""
    sampler = ModelSampler(model)
    for sensor in sensors.values():
        sampler._reference(sensor)
    cutoff = utc_time(as_of) if as_of is not None else None
    matches, excluded = [], []
    for obs in observations:
        sensor = sensors[obs.sensor_id]
        reason = ("missing_value" if obs.value is None else
                  "quality_flag" if obs.quality_flag not in include_quality else
                  "not_available_as_of" if cutoff is not None and obs.arrival_time > cutoff else None)
        if reason:
            excluded.append({"sensor_id": obs.sensor_id, "valid_time_utc": iso_utc(obs.valid_time), "reason": reason, "quality_flag": obs.quality_flag})
            continue
        predicted = sampler.predict(sensor, obs)
        matches.append({
            "sensor_id": sensor.sensor_id, "kind": sensor.kind, "role": sensor.role,
            "units": obs.units, "valid_time_utc": iso_utc(obs.valid_time),
            "arrival_time_utc": iso_utc(obs.arrival_time),
            "interval_start_utc": iso_utc(obs.interval_start) if obs.interval_start else None,
            "latency_seconds": (obs.arrival_time - obs.valid_time).total_seconds(),
            "observed": obs.value, "predicted": predicted, "uncertainty": obs.uncertainty,
            "residual_model_minus_observed": predicted - float(obs.value),
            "quality_flag": obs.quality_flag,
        })
    groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    by_sensor: dict[str, list[dict[str, Any]]] = {}
    for row in matches:
        groups.setdefault((row["kind"], row["units"], row["role"]), []).append(row)
        by_sensor.setdefault(row["sensor_id"], []).append(row)
    summaries = [{"kind": kind, "units": units, "role": role, **metrics(rows)}
                 for (kind, units, role), rows in sorted(groups.items())]
    return {
        "method": "Forward-operator comparison only; no assimilation or parameter fitting",
        "scope": "Synthetic self-consistency exercise; not independent field validation" if all("synthetic" in s.provenance.lower() for s in sensors.values()) else "Sensor comparison; independent calibration/validation provenance must be established separately",
        "state_and_discharge_validation": [row for row in summaries if row["kind"] != "precipitation"],
        "forcing_validation": [row for row in summaries if row["kind"] == "precipitation"],
        "by_sensor": [{"sensor_id": key, "kind": rows[0]["kind"], "units": rows[0]["units"], "role": rows[0]["role"], **metrics(rows)} for key, rows in sorted(by_sensor.items())],
        "matched_count": len(matches), "excluded_count": len(excluded),
        "matches": matches, "excluded": excluded,
        "as_of_utc": iso_utc(cutoff) if cutoff else None,
        "uncertainty_note": "Normalization uses observation 1-sigma uncertainty only; it is not a forecast confidence interval, calibrated innovation chi-square, or evidence of independent errors.",
        "operator_notes": [
            "State snapshots use interior trilinear point interpolation or axis-aligned support overlap; no extrapolation beyond grid edges.",
            "Screen means assume one saturated aquifer and uniform screen conductance per unit length; saturation is checked at supporting cell centres, so subcell water-table crossings need refinement or a richer well model.",
            "water_level_depth is local surface elevation minus measured hydraulic head; a confined potentiometric level is not a water table.",
            "Rate records are preceding-interval means and are integrated over each observation's exact interval.",
            "The outlet quantity is identified by model metadata; a fixed-head boundary flux is not independently validated river routing.",
        ],
    }


def write_synthetic_demo(model: dict[str, Any], destination: str | Path) -> tuple[Path, Path]:
    """Write explicitly synthetic observations sampled from an actual model run.

    Deterministic small offsets exercise residual metrics. They are not real
    measurement errors, and the sampled run is not an independent truth model.
    """
    sampler = ModelSampler(model)
    x, y, z = sampler.axes
    active = sampler.grid["active"]
    reference = sampler.meta["coordinate_reference"]
    choices = []
    for j in range(len(y)):
        for i in range(len(x)):
            levels = [k for k in range(len(z)) if active[k][j][i]]
            if len(levels) >= 3:
                choices.append((j, i, levels))
    if len(choices) < 2:
        raise ObservationError("Synthetic example needs two active columns with at least three cells")
    selected = [choices[len(choices) // 3], choices[2 * len(choices) // 3]]
    sensor_rows: list[dict[str, Any]] = []
    for n, (j, i, levels) in enumerate(selected):
        sensor_rows.append({"sensor_id": f"SM_{n + 1}", "kind": "soil_moisture", "x_m": x[i], "y_m": y[j], "z_m": z[levels[-2]], "units": "%", "uncertainty": 1.5, "role": "train" if n == 0 else "validation", "support_x_m": sampler.spacing[0], "support_y_m": sampler.spacing[1], "support_z_m": sampler.spacing[2]})
    j, i, levels = selected[1]
    sensor_rows.append({"sensor_id": "PSI_1", "kind": "pore_pressure", "x_m": x[i], "y_m": y[j], "z_m": z[levels[-2]], "units": "kPa", "uncertainty": 0.2, "role": "validation"})
    saturated = next(((j, i, levels) for j, i, levels in choices
                      if all(frame["pressure_head_m"][k][j][i] is not None and frame["pressure_head_m"][k][j][i] >= 0
                             for frame in sampler.snapshots for k in levels[:2])), None)
    if saturated is None:
        raise ObservationError("Synthetic head-well example needs a screen saturated throughout the run")
    j, i, levels = saturated
    sensor_rows.append({"sensor_id": "HEAD_1", "kind": "hydraulic_head", "x_m": x[i], "y_m": y[j], "screen_bottom_m": z[levels[0]] - sampler.spacing[2] / 2, "screen_top_m": z[levels[1]] + sampler.spacing[2] / 2, "aquifer_id": sampler.meta["aquifer_id"], "screen_assumption": "single_aquifer_uniform_conductance", "units": "m", "uncertainty": 0.03, "role": "validation"})
    outlet = sampler.meta["outlet"]
    sensor_rows.extend([
        {"sensor_id": "Q_OUT", "kind": "streamflow", "x_m": outlet.get("x_m", x[-1]), "y_m": outlet.get("y_m", y[-1]), "outlet_id": outlet["id"], "units": "L/s", "uncertainty": 0.02, "role": "validation"},
        {"sensor_id": "P_GAUGE", "kind": "precipitation", "x_m": x[len(x) // 2], "y_m": y[len(y) // 2], "units": "mm/h", "uncertainty": 0.03, "role": "validation"},
    ])
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    sensor_path, observation_path = destination / "sensors.csv", destination / "observations.csv"
    with sensor_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=SENSOR_COLUMNS)
        writer.writeheader()
        for row in sensor_rows:
            writer.writerow({**row, **reference, "provenance": "synthetic sampled from the same valley3d run; deterministic offsets; not field data"})
    sensors = load_sensors(sensor_path)
    end = min(sampler.times[-1], sampler.record_times[-1])
    sample_hours = sorted({round(end * fraction, 8) for fraction in (0.25, 0.5, 0.75, 1.0)})
    with observation_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OBSERVATION_COLUMNS)
        writer.writeheader()
        for n, hour in enumerate(sample_hours):
            valid = sampler.start + timedelta(hours=hour)
            for sensor in sensors.values():
                start = valid - timedelta(hours=min(1.0, hour)) if sensor.kind in RATE_KINDS else None
                obs = Observation(sensor.sensor_id, valid, valid + timedelta(minutes=10), None, sensor.uncertainty, sensor.units, "good", start)
                truth = sampler.predict(sensor, obs)
                factor = unit_factor(sensor.kind, sensor.units, (valid - start).total_seconds() if start else None)
                offset = sensor.uncertainty * (0.4 if n % 2 == 0 else -0.4)
                value = truth / factor + offset
                if sensor.kind in {"soil_moisture", "precipitation"}:
                    value = max(0.0, value)
                missing = sensor.sensor_id == "SM_1" and n == len(sample_hours) - 1
                writer.writerow({"sensor_id": sensor.sensor_id, "valid_time_utc": iso_utc(valid), "arrival_time_utc": iso_utc(valid + timedelta(minutes=10)), "interval_start_utc": iso_utc(start) if start else "", "value": "" if missing else f"{value:.10g}", "units": sensor.units, "uncertainty": sensor.uncertainty, "quality_flag": "missing" if missing else "good"})
    return sensor_path, observation_path


def main(argv: list[str] | None = None) -> int:
    """Command-line CSV ingestion and forward-comparison entry point."""
    parser = argparse.ArgumentParser(description="Ingest sensors and compare with precomputed model outputs; no assimilation")
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--sensors", type=Path)
    parser.add_argument("--observations", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--as-of", help="Only use data arrived by this timezone-aware time")
    parser.add_argument("--make-demo", "--generate-demo", dest="make_demo", type=Path,
                        help="Write synthetic example CSVs sampled from this model run to this directory")
    args = parser.parse_args(argv)
    try:
        model = json.loads(args.model.read_text())
        if args.make_demo is not None:
            args.sensors, args.observations = write_synthetic_demo(model, args.make_demo)
        if args.sensors is None or args.observations is None:
            parser.error("--sensors and --observations are required unless --make-demo is used")
        sensors = load_sensors(args.sensors)
        observations = load_observations(args.observations, sensors)
        report = compare_observations(model, sensors, observations, as_of=args.as_of)
        report["inputs"] = {"model": str(args.model), "sensors": str(args.sensors), "observations": str(args.observations)}
        report["input_sha256"] = {key: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                                  for key, path in report["inputs"].items()}
        report["model_coordinate_reference"] = model["meta"].get("coordinate_reference")
        report["model_outlet_definition"] = model["meta"].get("outlet")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(json.dumps({"output": str(args.output), "matched": report["matched_count"],
                          "excluded": report["excluded_count"], "method": report["method"]}))
        return 0
    except (ObservationError, OSError, ValueError, KeyError) as exc:
        print(f"Sensor validation failed: {exc}", file=sys.stderr)
        return 2
