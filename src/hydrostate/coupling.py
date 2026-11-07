"""Conservative fixed-grid material updates and detached hydraulic restarts.

The exchanged conserved quantity is the backend's storage variable
``W = theta(psi) + Ss max(psi, 0)``, not physical water content ``theta``.
These operations do not model deformation, air pressure, solid mass, damage
generation, heat, or seismic excess pore-pressure generation. They reconcile
prescribed material changes on an unchanged bulk-volume Cartesian domain.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import re
from typing import Any, Mapping

import numpy as np

from hydrostate.flow import Grid3D, RichardsModel


PROTOCOL = "hydrostate.hydraulic-checkpoint/1"
DEFAULT_PRESSURE_BOUNDS_M = (-100_000.0, 100.0)
MATERIAL_NAMES = ("theta_r", "theta_s", "alpha_per_m", "n", "ks_m_s", "ss_per_m")
MATERIAL_UNITS = {"theta_r": "m3/m3", "theta_s": "m3/m3", "alpha_per_m": "1/m",
                  "n": "1", "ks_m_s": "m/s", "ss_per_m": "1/m", "linear_capacity": "1/m"}
STORAGE_TOLERANCE = 5e-13


class CouplingError(ValueError):
    """A proposed interchange cannot preserve the stated hydraulic contract."""


def _utc(value: str | datetime, name: str) -> str:
    if isinstance(value, str) and re.search(r"[.,]\d{7,}", value):
        raise CouplingError(f"{name}: hydraulic clock supports microsecond resolution; sub-microsecond precision is rejected")
    try:
        time = datetime.fromisoformat(value.replace("Z", "+00:00")) if isinstance(value, str) else value
        if not isinstance(time, datetime) or time.tzinfo is None or time.utcoffset() is None:
            raise ValueError("Timezone required")
        if getattr(time, "nanosecond", 0) != 0:
            raise CouplingError(f"{name}: hydraulic clock supports microsecond resolution; nonzero nanoseconds are rejected")
        normalized = time.astimezone(timezone.utc)
        # Strip datetime subclasses only after checking that no precision is lost.
        normalized = datetime(normalized.year, normalized.month, normalized.day, normalized.hour,
                              normalized.minute, normalized.second, normalized.microsecond, tzinfo=timezone.utc)
        return normalized.isoformat().replace("+00:00", "Z")
    except CouplingError:
        raise
    except (ValueError, TypeError) as exc:
        raise CouplingError(f"{name} requires a timezone-aware ISO-8601 timestamp") from exc


def _reference(value: Mapping[str, Any]) -> dict:
    if not isinstance(value, Mapping):
        raise CouplingError("reference requires an object")
    result = deepcopy(dict(value))
    for name in ("horizontal_crs", "vertical_datum"):
        if not isinstance(result.get(name), str) or not result[name].strip():
            raise CouplingError(f"reference.{name} requires a nonempty string")
    if result.get("length_unit") != "m" or result.get("z_positive") != "up":
        raise CouplingError("reference requires length_unit='m' and z_positive='up'")
    return result


def _json_copy(value: Any) -> Any:
    try:
        return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))
    except (TypeError, ValueError) as exc:
        raise CouplingError("Checkpoint metadata must contain standard finite JSON values") from exc


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _field(value: Any, grid: Grid3D, name: str) -> np.ndarray:
    try:
        raw = np.asarray(value, dtype=object)
        if raw.shape != grid.shape:
            raise CouplingError(f"{name} must have shape {grid.shape} [z,y,x]")
        if any(v is None or isinstance(v, (str, bool)) for v in raw[grid.active]):
            raise CouplingError(f"{name} active cells require finite numbers")
        result = raw.astype(float)
    except (TypeError, ValueError) as exc:
        raise CouplingError(f"{name} requires finite numeric active values") from exc
    if not np.all(np.isfinite(result[grid.active])):
        raise CouplingError(f"{name} active cells require finite numbers")
    inactive = result[~grid.active]
    if np.any(~np.isnan(inactive)):
        raise CouplingError(f"{name} inactive cells must be null/NaN, not assigned material or water")
    return result.copy()


def _json_field(value: np.ndarray) -> list:
    return np.where(np.isfinite(value), value, None).tolist()


def grid_payload(grid: Grid3D) -> dict:
    """Detached SI Cartesian geometry, preserving the exact active mask."""
    return {"x_m": grid.x_m.tolist(), "y_m": grid.y_m.tolist(), "z_m": grid.z_m.tolist(),
            "dx_m": float(grid.dx_m), "dy_m": float(grid.dy_m), "dz_m": float(grid.dz_m),
            "active": grid.active.tolist()}


def grid_fingerprint(grid: Grid3D) -> str:
    return _digest(grid_payload(grid))


def _expanded(model: RichardsModel, values: np.ndarray) -> np.ndarray:
    result = np.full(model.grid.shape, np.nan)
    result.ravel()[model.ids] = values
    return result


def material_payload(model: RichardsModel) -> dict:
    """Detached complete backend material state, including verification mode."""
    output = {name: {"units": MATERIAL_UNITS[name], "values": _json_field(_expanded(model, getattr(model, name)))}
              for name in MATERIAL_NAMES}
    output["linear_capacity"] = (None if model.linear_capacity is None else
                                 {"units": "1/m", "values": _json_field(_expanded(model, model.linear_capacity))})
    return output


def model_fingerprint(model: RichardsModel) -> str:
    """Hash actual geometry/material values; detects in-place material mutation."""
    return _digest({"grid": grid_payload(model.grid), "materials": material_payload(model)})


def _bounds(value: tuple[float, float]) -> tuple[float, float]:
    try:
        if len(value) != 2 or any(isinstance(v, bool) for v in value):
            raise ValueError("Two finite bounds required")
        lo, hi = map(float, value)
    except (TypeError, ValueError) as exc:
        raise CouplingError("pressure_bounds_m must contain two finite increasing bounds") from exc
    if not np.isfinite([lo, hi]).all() or lo >= hi:
        raise CouplingError("pressure_bounds_m must contain two finite increasing bounds")
    return lo, hi


def _check_pressure(pressure: np.ndarray, grid: Grid3D, bounds: tuple[float, float]) -> None:
    lo, hi = _bounds(bounds)
    if np.any((pressure[grid.active] < lo) | (pressure[grid.active] > hi)):
        raise CouplingError(f"Pressure head exceeds declared admissible envelope [{lo}, {hi}] m")


def pressure_from_storage(
    model: RichardsModel, stored_water: Any, *, saturated_pressure_head_m: Any = None,
    pressure_bounds_m: tuple[float, float] = DEFAULT_PRESSURE_BOUNDS_M,
) -> np.ndarray:
    """Invert actual model storage without clipping, extrapolation or lost water.

    For ``Ss=0, W=theta_s`` saturated pressure is nonunique: supply an explicit
    nonnegative head hint. The default pressure envelope is a synthetic guard,
    not a calibrated hydrogeologic or mechanical operating range.
    """
    water = _field(stored_water, model.grid, "stored_water")
    w = water.ravel()[model.ids]
    if model.linear_capacity is not None:
        pressure = w / model.linear_capacity
    else:
        tr, ts, ss = model.theta_r, model.theta_s, model.ss_per_m
        if np.any(w <= tr):
            raise CouplingError("W <= theta_r has no finite equilibrium pressure representation")
        pressure = np.empty_like(w)
        unsaturated = w < ts
        se = (w[unsaturated] - tr[unsaturated]) / (ts[unsaturated] - tr[unsaturated])
        m = 1 - 1 / model.n[unsaturated]
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            pressure[unsaturated] = -np.expm1(-np.log(se) / m) ** (1 / model.n[unsaturated]) / model.alpha_per_m[unsaturated]
        compressible = ~unsaturated & (ss > 0)
        pressure[compressible] = (w[compressible] - ts[compressible]) / ss[compressible]
        nonunique = ~unsaturated & (ss == 0)
        if np.any(w[nonunique] > ts[nonunique]):
            raise CouplingError("W > theta_s with Ss=0 exceeds storage capacity")
        if np.any(nonunique):
            if saturated_pressure_head_m is None:
                raise CouplingError("Ss=0 saturated storage cannot determine pressure; supply a nonnegative pressure hint")
            hint = _field(saturated_pressure_head_m, model.grid, "saturated_pressure_head_m").ravel()[model.ids]
            if np.any(hint[nonunique] < 0):
                raise CouplingError("Ss=0 saturated pressure hint must be nonnegative")
            pressure[nonunique] = hint[nonunique]
    if not np.all(np.isfinite(pressure)):
        raise CouplingError("Stored water requires nonfinite pressure under the supplied retention law")
    full = _expanded(model, pressure)
    _check_pressure(full, model.grid, pressure_bounds_m)
    error = model.storage(full)[model.grid.active] - water[model.grid.active]
    if np.max(abs(error), initial=0) > STORAGE_TOLERANCE:
        raise CouplingError("Storage inverse failed its per-cell round-trip tolerance")
    return full


@dataclass(frozen=True)
class MaterialUpdateResult:
    pressure_head_m: np.ndarray
    storage_before: np.ndarray
    storage_after: np.ndarray
    delta_theta: np.ndarray
    delta_storage: np.ndarray
    signed_residual_m3: float
    absolute_residual_m3: float
    max_abs_cell_residual_m3: float


def reconcile_material_change(
    old_model: RichardsModel, pressure_head_m: Any, new_model: RichardsModel, *,
    stored_water: Any = None, saturated_pressure_head_m: Any = None,
    pressure_bounds_m: tuple[float, float] = DEFAULT_PRESSURE_BOUNDS_M,
) -> MaterialUpdateResult:
    """Preserve each cell's W when replacing fixed-grid material parameters.

    No object is mutated. Build ``new_model`` separately and install it only
    after this function succeeds. A checkpoint's explicit ``stored_water``
    detects a changed old material law instead of silently recomputing mass.
    The unchanged old pressure is the declared default hint for an ambiguous
    new ``Ss=0`` saturated branch, and is accepted only if nonnegative.
    """
    if grid_fingerprint(old_model.grid) != grid_fingerprint(new_model.grid):
        raise CouplingError("Grid/active-mask change requires explicit conservative remeshing and geometry mechanics")
    if (old_model.linear_capacity is None) != (new_model.linear_capacity is None):
        raise CouplingError("Cannot exchange physical water and linear-verification storage modes")
    psi = _field(pressure_head_m, old_model.grid, "pressure_head_m")
    _check_pressure(psi, old_model.grid, pressure_bounds_m)
    computed = old_model.storage(psi)
    water = computed.copy() if stored_water is None else _field(stored_water, old_model.grid, "stored_water")
    if np.max(abs((computed - water)[old_model.grid.active]), initial=0) > STORAGE_TOLERANCE:
        raise CouplingError("Old pressure/material law disagrees with explicit stored water; possible material mutation")
    hint = psi if saturated_pressure_head_m is None else saturated_pressure_head_m
    # Preserve exactly equal pressure when the entire storage law is unchanged;
    # Ks-only updates should not incur inversion roundoff or lose Ss=0 pressure.
    storage_names = ("theta_r", "theta_s", "alpha_per_m", "n", "ss_per_m", "linear_capacity")
    same_storage = all(np.array_equal(getattr(old_model, name), getattr(new_model, name)) for name in storage_names)
    updated = psi.copy() if same_storage else pressure_from_storage(
        new_model, water, saturated_pressure_head_m=hint, pressure_bounds_m=pressure_bounds_m)
    after = new_model.storage(updated)
    difference = after - water
    local_m3 = difference[old_model.grid.active] * old_model.grid.cell_volume_m3
    return MaterialUpdateResult(updated, water.copy(), after, new_model.theta(updated) - old_model.theta(psi),
                                difference, float(local_m3.sum()), float(abs(local_m3).sum()),
                                float(np.max(abs(local_m3), initial=0)))


def build_checkpoint(
    model: RichardsModel, pressure_head_m: Any, *, valid_time_utc: str | datetime,
    reference: Mapping[str, Any], issued_time_utc: str | datetime | None = None,
    metadata: Mapping[str, Any] | None = None,
    pressure_bounds_m: tuple[float, float] = DEFAULT_PRESSURE_BOUNDS_M,
) -> dict:
    """Detach complete state of this equilibrium Richards component only.

    Future snow, channel, hysteresis, mechanics/event memory and forcing state
    need their own owners' checkpoints in an enclosing coupled restart.
    """
    psi = _field(pressure_head_m, model.grid, "pressure_head_m")
    _check_pressure(psi, model.grid, pressure_bounds_m)
    result = {"protocol": PROTOCOL, "array_order": "zyx", "reference": _reference(reference),
              "valid_time_utc": _utc(valid_time_utc, "valid_time_utc"),
              "issued_time_utc": _utc(issued_time_utc or datetime.now(timezone.utc), "issued_time_utc"),
              "grid": grid_payload(model.grid), "materials": material_payload(model),
              "state": {"pressure_head_m": {"units": "m", "values": _json_field(psi)},
                        "stored_water": {"units": "m3/m3", "definition": "W=theta+Ss*max(psi,0); C*psi in linear verification mode",
                                         "values": _json_field(model.storage(psi))},
                        "theta": {"units": "m3/m3", "values": _json_field(model.theta(psi))}},
              "pressure_bounds_m": list(_bounds(pressure_bounds_m)),
              "grid_sha256": grid_fingerprint(model.grid), "model_sha256": model_fingerprint(model),
              "metadata": _json_copy(dict(metadata or {}))}
    result["content_sha256"] = _digest(result)
    return _json_copy(result)


@dataclass(frozen=True)
class RestoredCheckpoint:
    model: RichardsModel
    pressure_head_m: np.ndarray
    stored_water: np.ndarray
    valid_time_utc: str
    reference: dict
    metadata: dict


def _restore_grid(value: Any) -> Grid3D:
    if not isinstance(value, dict):
        raise CouplingError("Checkpoint grid requires an object")
    try:
        raw = np.asarray(value["active"], dtype=object)
        if any(type(v) is not bool for v in raw.ravel()):
            raise CouplingError("Checkpoint active mask requires JSON booleans")
        for name in ("x_m", "y_m", "z_m"):
            if any(isinstance(v, (str, bool)) or v is None for v in np.asarray(value[name], dtype=object).ravel()):
                raise CouplingError("Checkpoint coordinates require numeric SI values")
        for name in ("dx_m", "dy_m", "dz_m"):
            if isinstance(value[name], (str, bool)) or not isinstance(value[name], (int, float)):
                raise CouplingError("Checkpoint spacings require finite positive numeric SI values")
        return Grid3D(value["x_m"], value["y_m"], value["z_m"],
                      value["dx_m"], value["dy_m"], value["dz_m"], raw.astype(bool))
    except (KeyError, TypeError, ValueError) as exc:
        raise CouplingError(f"Invalid checkpoint grid: {exc}") from exc


def restore_checkpoint(
    payload: Mapping[str, Any], *, expected_reference: Mapping[str, Any] | None = None,
    expected_grid: Grid3D | None = None, expected_model: RichardsModel | None = None,
    expected_valid_time_utc: str | datetime | None = None,
) -> RestoredCheckpoint:
    """Validate integrity and rebuild a detached backend with its actual water.

    Optional expected values are mandatory when restoring into an existing
    known domain. No resampling, unit conversion or clock interpolation occurs.
    A content hash is an integrity check, not a cryptographic authenticity claim.
    """
    p = _json_copy(payload)
    if not isinstance(p, dict) or p.get("protocol") != PROTOCOL or p.get("array_order") != "zyx":
        raise CouplingError(f"Checkpoint requires protocol={PROTOCOL!r}, array_order='zyx'")
    content_hash = p.pop("content_sha256", None)
    if content_hash != _digest(p):
        raise CouplingError("Checkpoint content hash mismatch")
    reference = _reference(p.get("reference"))
    if expected_reference is not None and reference != _reference(expected_reference):
        raise CouplingError("Checkpoint coordinate reference/datum mismatch")
    valid = _utc(p.get("valid_time_utc"), "valid_time_utc")
    _utc(p.get("issued_time_utc"), "issued_time_utc")
    if expected_valid_time_utc is not None and valid != _utc(expected_valid_time_utc, "expected_valid_time_utc"):
        raise CouplingError("Checkpoint valid-time clock mismatch")
    grid = _restore_grid(p.get("grid"))
    if grid_fingerprint(grid) != p.get("grid_sha256"):
        raise CouplingError("Checkpoint grid fingerprint mismatch")
    if expected_grid is not None and grid_fingerprint(grid) != grid_fingerprint(expected_grid):
        raise CouplingError("Checkpoint grid/active-mask mismatch")
    materials = p.get("materials")
    if not isinstance(materials, dict) or set(materials) != set(MATERIAL_NAMES) | {"linear_capacity"}:
        raise CouplingError("Checkpoint requires complete backend materials")
    arrays = {}
    for name, descriptor in materials.items():
        if name == "linear_capacity" and descriptor is None:
            arrays[name] = None
            continue
        if not isinstance(descriptor, dict) or descriptor.get("units") != MATERIAL_UNITS[name]:
            raise CouplingError(f"Checkpoint material {name} units mismatch")
        arrays[name] = _field(descriptor.get("values"), grid, name)
    try:
        model = RichardsModel(grid, **arrays)
    except ValueError as exc:
        raise CouplingError(f"Invalid checkpoint hydraulic materials: {exc}") from exc
    if model_fingerprint(model) != p.get("model_sha256"):
        raise CouplingError("Checkpoint model/material fingerprint mismatch")
    if expected_model is not None and model_fingerprint(model) != model_fingerprint(expected_model):
        raise CouplingError("Live model/material fingerprint mismatch; a material may have been mutated")
    fields = p.get("state")
    if not isinstance(fields, dict) or set(fields) != {"pressure_head_m", "stored_water", "theta"}:
        raise CouplingError("Checkpoint requires pressure, physical theta and conserved W")
    state = {}
    for name, units in (("pressure_head_m", "m"), ("stored_water", "m3/m3"), ("theta", "m3/m3")):
        descriptor = fields[name]
        if not isinstance(descriptor, dict) or descriptor.get("units") != units:
            raise CouplingError(f"Checkpoint state {name} units mismatch")
        state[name] = _field(descriptor.get("values"), grid, name)
    _check_pressure(state["pressure_head_m"], grid, p.get("pressure_bounds_m"))
    for name, calculated in (("stored_water", model.storage(state["pressure_head_m"])),
                             ("theta", model.theta(state["pressure_head_m"]))):
        if np.max(abs((state[name] - calculated)[grid.active]), initial=0) > STORAGE_TOLERANCE:
            raise CouplingError(f"Checkpoint {name} disagrees with pressure/material constitutive law")
    return RestoredCheckpoint(model, state["pressure_head_m"].copy(), state["stored_water"].copy(),
                              valid, reference, deepcopy(p.get("metadata", {})))
