"""Typed model interchange through Gaia's actual xarray/DataTree Zarr writer.

This is a model-output contract layered on Gaia's storage contract. Gaia's
horizontal raster alignment and temporal averaging are not model coupling or
restart operations. Synthetic local coordinates deliberately carry no CRS.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import hashlib
import json
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import xarray as xr


GAIA_COMMIT = "aa8ed3def5af4c32a4a1e601252fcee439344ac8"
CONTRACT = "hydrostate.mechanicalstate.datahub/1"
STATE_GROUP = "rasters/hydrostate-state"
STATIC_GROUP = "rasters/mechanicalstate-materials"
MATERIAL_GROUP = "rasters/mechanicalstate-state"
INTERVAL_GROUP = "timeseries/hydrostate-fluxes"
EVENT_GROUP = "events/material-jumps"
CHECKPOINT_GROUP = "checkpoints/hydro-mechanical"
SPATIAL_DIMS = ("z", "y", "x")

# This is an interface registry, not a CF standard-name table. No standard_name
# is fabricated for suction, seismic damage, water-table proxies or materials.
FIELD_UNITS = {
    "pressure_head_m": "m", "pressure_head_before_m": "m", "pressure_head_after_m": "m",
    "hydraulic_head_m": "m", "total_head_m": "m", "pore_pressure_pa": "Pa",
    "theta": "m3/m3", "water_content": "m3/m3", "theta_r": "m3/m3", "theta_s": "m3/m3",
    "water_storage_fraction": "m3/m3", "storage_water_equivalent": "m3/m3",
    "total_water_storage_fraction": "m3/m3", "stored_water": "m3/m3",
    "porosity": "1", "phi": "1", "saturation": "1", "water_saturation": "1",
    "damage": "1", "burn_fraction": "1", "water_repellency": "1", "root_loss": "1",
    "elastic_damage": "1", "logK": "1", "logalpha": "1",
    "log_conductivity_anomaly": "1", "log_alpha_anomaly": "1",
    "water_storage_before": "m3/m3", "water_storage_after": "m3/m3",
    "alpha_per_m": "1/m", "n": "1", "ss_per_m": "1/m", "ks_m_s": "m/s",
    "intrinsic_permeability_m2": "m2", "vp_m_s": "m/s", "vs_m_s": "m/s",
    "bulk_density_kg_m3": "kg/m3", "facies_id": "1", "active": "1",
}


@dataclass(frozen=True)
class Field:
    """Array with explicit units and optional quantity/support metadata.

    Floating-point width, integer width and booleans are retained. Object
    arrays are accepted only when every value is a string; missing strings
    must be represented by an explicit quality flag rather than None.
    """

    values: Any
    units: str
    attrs: Mapping[str, Any] = field(default_factory=dict)


def _jsonable(value: Any) -> Any:
    try:
        return json.loads(json.dumps(value, allow_nan=False))
    except (TypeError, ValueError) as exc:
        raise ValueError("Metadata must be finite JSON-serializable values") from exc


def _utc(values: Sequence[Any], name: str, *, increasing: bool = True) -> np.ndarray:
    result = []
    for value in values:
        if not isinstance(value, (datetime, str)):
            raise ValueError(f"{name} requires timezone-aware datetime or ISO-8601 text")
        try:
            parsed = pd.Timestamp(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} requires timezone-aware ISO-8601 timestamps") from exc
        if pd.isna(parsed) or parsed.tzinfo is None or parsed.utcoffset() is None:
            raise ValueError(f"{name} requires a timezone, preferably UTC")
        # pandas preserves nanoseconds in ISO strings, unlike datetime.fromisoformat.
        result.append(parsed.tz_convert("UTC").tz_localize(None).as_unit("ns").to_datetime64())
    array = np.asarray(result, dtype="datetime64[ns]")
    if not array.size or np.isnat(array).any():
        raise ValueError(f"{name} cannot be empty or contain NaT")
    if increasing and array.size > 1 and np.any(np.diff(array) <= np.timedelta64(0, "ns")):
        raise ValueError(f"{name} must be unique and increasing")
    return array


def _iso(value: np.datetime64) -> str:
    return np.datetime_as_string(value, unit="ns") + "Z"


def source_attrs(*, title: str, description: str, source_commit: str,
                 lineage: str, repository: str = "https://github.com/gaia-hazlab/hydrostate",
                 license: str = "MIT") -> dict:
    """Minimal Gaia source metadata for synthetic/model output, not field data."""
    return {
        "title": title, "description": description, "license": license,
        "providers": [{"name": "GAIA Hazlab", "roles": ["producer", "processor"],
                       "url": "https://github.com/gaia-hazlab"}],
        "links": [{"rel": "derived_from", "href": f"{repository}/tree/{source_commit}"},
                  {"rel": "license", "href": f"{repository}/blob/{source_commit}/LICENSE"}],
        "processing:lineage": lineage, "version": source_commit,
    }


def _source_attrs(provenance: Mapping[str, Any]) -> dict:
    result = _jsonable(dict(provenance))
    for name in ("title", "description", "license", "processing:lineage"):
        if not isinstance(result.get(name), str) or not result[name].strip():
            raise ValueError(f"Source provenance requires nonempty {name}")
    for name in ("providers", "links"):
        if not isinstance(result.get(name), list) or not result[name]:
            raise ValueError(f"Source provenance requires {name}")
    if not any(link.get("rel") in ("derived_from", "via") for link in result["links"]):
        raise ValueError("Source links must identify a derived_from or via source")
    if not any(link.get("rel") == "license" for link in result["links"]):
        raise ValueError("Source links must identify license terms")
    return result


def _reference(reference: Mapping[str, Any]) -> tuple[dict, Any]:
    result = _jsonable(dict(reference))
    if result.get("length_unit") != "m" or result.get("z_positive") != "up":
        raise ValueError("Reference requires length_unit='m' and z_positive='up'")
    datum = result.get("vertical_datum")
    horizontal = result.get("horizontal_crs")
    if not isinstance(datum, str) or not datum.strip() or not isinstance(horizontal, str) or not horizontal.strip():
        raise ValueError("Horizontal CRS/local reference and vertical datum are required")
    local = horizontal.upper().startswith("LOCAL")
    if "synthetic_local" in result and type(result["synthetic_local"]) is not bool:
        raise ValueError("synthetic_local must be an actual boolean")
    if result.get("synthetic_local", local) != local:
        raise ValueError("Synthetic-local status and horizontal_crs disagree")
    result["synthetic_local"] = local
    if local:
        return result, None
    import pyproj
    crs = pyproj.CRS.from_user_input(horizontal)
    if not crs.is_projected or any(abs(axis.unit_conversion_factor - 1) > 1e-12 for axis in crs.axis_info[:2]):
        raise ValueError("Cartesian solver x/y require a projected metre CRS; transform explicitly upstream")
    return result, crs


def canonical_reference(reference: Mapping[str, Any]) -> dict:
    """Validate and normalize the shared metric reference without transforming it."""
    normalized, _ = _reference(reference)
    return normalized


def _hash(parts: Sequence[Any]) -> str:
    digest = hashlib.sha256()
    for value in parts:
        if isinstance(value, np.ndarray):
            digest.update(("utf8" if value.dtype.kind in "OUST" else str(value.dtype)).encode())
            digest.update(str(value.shape).encode())
            if value.dtype.kind in "OUST":
                digest.update(json.dumps(value.tolist(), ensure_ascii=False, separators=(",", ":")).encode())
            else:
                digest.update(np.ascontiguousarray(value).tobytes())
        else:
            digest.update(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())
    return digest.hexdigest()


def _array(name: str, descriptor: Field, shape: tuple[int, ...]) -> tuple[np.ndarray, dict]:
    if not isinstance(descriptor, Field) or not isinstance(descriptor.units, str) or not descriptor.units:
        raise ValueError(f"{name} requires a Field with explicit units")
    if name in FIELD_UNITS and descriptor.units != FIELD_UNITS[name]:
        raise ValueError(f"{name} requires units={FIELD_UNITS[name]!r}; convert explicitly")
    array = np.asarray(descriptor.values)
    if array.shape != shape:
        raise ValueError(f"{name} requires shape {shape}, received {array.shape}")
    if array.dtype.kind not in "fiuUbOST":
        raise ValueError(f"{name} has unsupported dtype {array.dtype}")
    if array.dtype.kind == "O" and not all(isinstance(value, str) for value in array.flat):
        raise ValueError(f"{name}: object arrays must contain only text")
    if array.dtype.kind == "f" and np.isinf(array).any():
        raise ValueError(f"{name} contains infinity; use NaN for missing data")
    attrs = _jsonable(dict(descriptor.attrs))
    if "standard_name" in attrs:
        raise ValueError("Custom field standard_name is not accepted; use long_name and units")
    if "units" in attrs and attrs["units"] != descriptor.units:
        raise ValueError(f"{name} unit metadata disagrees")
    attrs.update(units=descriptor.units)
    return np.array(array, copy=True), attrs


def _data_vars(fields: Mapping[str, Field], dims: tuple[str, ...], shape: tuple[int, ...]) -> dict:
    result = {}
    for name, descriptor in fields.items():
        if not isinstance(name, str) or not name or name in (*SPATIAL_DIMS, "time", "event", "spatial_ref", "time_bounds", "active"):
            raise ValueError(f"Reserved/invalid field name {name!r}")
        array, attrs = _array(name, descriptor, shape)
        result[name] = (dims, array, attrs)
    return result


def _physical_domain(name: str, values: np.ndarray, active: np.ndarray) -> None:
    """A prognostic/material field is complete in soil and absent in air."""
    if name not in FIELD_UNITS or name in ("facies_id", "active") or values.ndim not in (3, 4):
        return
    if values.dtype.kind != "f":
        raise ValueError(f"Physical continuous field {name} requires floating values")
    mask = active if values.ndim == 3 else np.broadcast_to(active, values.shape)
    if not np.isfinite(values[mask]).all():
        raise ValueError(f"{name} requires finite values in every active cell")
    if not np.isnan(values[~mask]).all():
        raise ValueError(f"{name} requires NaN outside the active material domain")


def build_tree(*, grid: Any, reference: Mapping[str, Any], provenance: Mapping[str, Any],
               state_times_utc: Sequence[Any], state_fields: Mapping[str, Field],
               static_fields: Mapping[str, Field] | None = None,
               material_fields: Mapping[str, Field] | None = None,
               material_epochs: Sequence[Any] | None = None,
               interval_bounds_utc: Sequence[Sequence[Any]] | None = None,
               interval_fields: Mapping[str, Field] | None = None,
               event_times_utc: Sequence[Any] | None = None,
               event_ids: Sequence[str] | None = None,
               event_fields: Mapping[str, Field] | None = None,
               checkpoint_payload: Mapping[str, Any] | None = None) -> xr.DataTree:
    """Build exact-grid model groups without regridding or temporal averaging.

    State/dynamic material fields are ``(time,z,y,x)``; static fields ``(z,y,x)``.
    Interval/event fields may be a scalar ledger per interval/event or a whole
    ``(z,y,x)`` field per interval/event. Event left/right states are separate
    explicitly named variables; event times may repeat only with unique IDs.
    Naive datetimes and datetime64 input are rejected because neither declares
    a timezone. The stored CF datetime64 coordinates explicitly declare UTC.
    """
    if not state_fields:
        raise ValueError("At least one model state field is required")
    ref, crs = _reference(reference)
    source = _source_attrs(provenance)
    times = _utc(state_times_utc, "state_times_utc")
    xyz = {axis: np.asarray(getattr(grid, f"{axis}_m")) for axis in "xyz"}
    shape = tuple(len(xyz[axis]) for axis in "zyx")
    active = np.asarray(grid.active)
    if active.shape != shape or active.dtype.kind != "b" or not active.any():
        raise ValueError("Grid requires a boolean active mask matching (z,y,x)")
    for axis in "xyz":
        coord = xyz[axis]
        spacing = getattr(grid, f"d{axis}_m")
        if (coord.ndim != 1 or not coord.size or not np.isfinite(coord).all()
                or not np.isfinite(spacing) or spacing <= 0
                or (coord.size > 1 and not np.allclose(np.diff(coord), spacing, rtol=1e-10, atol=1e-12))):
            raise ValueError(f"Grid {axis} requires ascending uniform finite cell centers in metres")
    grid_hash = _hash([xyz["x"], xyz["y"], xyz["z"], active,
                       {axis: float(getattr(grid, f"d{axis}_m")) for axis in "xyz"}])
    reference_hash = _hash([ref])
    epochs = np.asarray(["static"] * len(times) if material_epochs is None else material_epochs)
    if epochs.shape != (len(times),) or epochs.dtype.kind not in "iuUOST":
        raise ValueError("material_epochs must have one integer/text value per state time")
    if epochs.dtype.kind == "O" and not all(isinstance(value, str) for value in epochs):
        raise ValueError("Object material_epochs must contain only text")
    coords = {axis: (axis, values.copy(), {"units": "m", "axis": axis.upper(),
                                          **({"positive": "up", "vertical_datum": ref["vertical_datum"]} if axis == "z" else {})})
              for axis, values in xyz.items()}
    coords["time"] = ("time", times, {"standard_name": "time", "axis": "T", "timezone": "UTC"})
    attrs = {**source, "hydrostate:contract": CONTRACT, "hydrostate:reference": ref,
             "hydrostate:grid-fingerprint": grid_hash, "hydrostate:reference-fingerprint": reference_hash,
             "hydrostate:gaia-contract-commit": GAIA_COMMIT,
             "hydrostate:alignment-policy": "no automatic gaia align on model states/materials/events",
             "start_datetime": _iso(times[0]), "end_datetime": _iso(times[-1])}

    def spatial(ds: xr.Dataset) -> xr.Dataset:
        ds["active"] = (SPATIAL_DIMS, active.copy(), {"units": "1", "long_name": "active physical storage cell"})
        ds["active"].encoding["_FillValue"] = None
        for axis in "xyz":
            ds.attrs[f"hydrostate:d{axis}_m"] = float(getattr(grid, f"d{axis}_m"))
        if crs is not None:
            import rioxarray  # noqa: F401: registers rio
            ds = ds.rio.write_crs(crs)
            ds = ds.rio.write_transform()
        else:
            ds.attrs["hydrostate:spatial-status"] = "synthetic local Cartesian; no earth georeferencing"
        return ds

    states = xr.Dataset(_data_vars(state_fields, ("time", *SPATIAL_DIMS), (len(times), *shape)),
                        coords=coords, attrs={**attrs, "hydrostate:role": "instantaneous prognostic state"})
    states = spatial(states).assign_coords(material_epoch=("time", epochs.copy()))
    for name in state_fields:
        _physical_domain(name, states[name].values, active)
        states[name].attrs.setdefault("cell_methods", "time: point")
    nodes = {"/": xr.Dataset(attrs={"title": source["title"], "hydrostate:contract": CONTRACT,
                                   "hydrostate:gaia-contract-commit": GAIA_COMMIT,
                                   "description": "Coupled model output; state and event groups require model-aware readers"}),
             STATE_GROUP: states}
    static = static_fields or {}
    dynamic = material_fields or {}
    static_arrays = _data_vars(static, SPATIAL_DIMS, shape)
    dynamic_arrays = _data_vars(dynamic, ("time", *SPATIAL_DIMS), (len(times), *shape))
    for name, (_, array, _) in {**static_arrays, **dynamic_arrays}.items():
        _physical_domain(name, array, active)
    fingerprints = []
    for i in range(len(times)):
        parts: list[Any] = []
        for name in sorted(static_arrays):
            _, array, field_attrs = static_arrays[name]
            parts.extend([name, field_attrs["units"], array])
        for name in sorted(dynamic_arrays):
            _, array, field_attrs = dynamic_arrays[name]
            parts.extend([name, field_attrs["units"], array[i]])
        fingerprints.append(_hash(parts))
    nodes[STATE_GROUP] = states.assign_coords(material_fingerprint=("time", np.asarray(fingerprints)))
    if static:
        nodes[STATIC_GROUP] = spatial(xr.Dataset(static_arrays,
            coords={axis: coords[axis] for axis in "xyz"},
            attrs={**attrs, "hydrostate:role": "reference material properties"}))
    if dynamic:
        materials = spatial(xr.Dataset(dynamic_arrays, coords=coords,
                                       attrs={**attrs, "hydrostate:role": "material state at matching valid times"}))
        nodes[MATERIAL_GROUP] = materials.assign_coords(material_epoch=("time", epochs.copy()),
                                                       material_fingerprint=("time", np.asarray(fingerprints)))
    if interval_fields or interval_bounds_utc is not None:
        if not interval_fields or interval_bounds_utc is None:
            raise ValueError("Interval fields require explicit bounds, and bounds require fields")
        bounds = np.stack([_utc(pair, "interval bounds") for pair in interval_bounds_utc])
        if bounds.ndim != 2 or bounds.shape[1] != 2 or np.any(bounds[:, 1] <= bounds[:, 0]):
            raise ValueError("Each interval requires increasing (start,end) timestamps")
        if len(bounds) > 1 and np.any(bounds[1:, 0] < bounds[:-1, 1]):
            raise ValueError("Intervals must be ordered and cannot overlap")
        values = {}
        for name, descriptor in interval_fields.items():
            is_grid = np.asarray(descriptor.values).ndim == 4
            dims = ("time", *SPATIAL_DIMS) if is_grid else ("time",)
            expected = (len(bounds), *shape) if is_grid else (len(bounds),)
            values.update(_data_vars({name: descriptor}, dims, expected))
            _physical_domain(name, values[name][1], active)
            method = values[name][2].get("cell_methods")
            if method not in ("time: mean", "time: sum"):
                raise ValueError(f"Interval field {name} requires cell_methods time: mean or time: sum")
        interval_coords = {"time": ("time", bounds[:, 1], {"standard_name": "time", "axis": "T", "timezone": "UTC", "bounds": "time_bounds"})}
        if any(len(value[0]) > 1 for value in values.values()):
            interval_coords.update({axis: coords[axis] for axis in "xyz"})
        intervals = xr.Dataset(values, coords=interval_coords,
            attrs={**attrs, "hydrostate:role": "integrated ledger / interval fluxes",
                   "start_datetime": _iso(bounds[0, 0]), "end_datetime": _iso(bounds[-1, 1])})
        intervals["time_bounds"] = (("time", "nv"), bounds)
        intervals["time_bounds"].attrs["timezone"] = "UTC"
        if "x" in intervals.dims:
            intervals = spatial(intervals)
        nodes[INTERVAL_GROUP] = intervals
    if event_fields or event_times_utc is not None or event_ids is not None:
        if not event_fields or event_times_utc is None or event_ids is None:
            raise ValueError("Events require exact valid times, unique IDs and fields")
        event_times = _utc(event_times_utc, "event_times_utc", increasing=False)
        if len(event_times) > 1 and np.any(np.diff(event_times) < np.timedelta64(0, "ns")):
            raise ValueError("Events must be ordered; simultaneous events retain their ID order")
        ids = np.asarray(event_ids)
        if ids.shape != event_times.shape or ids.dtype.kind not in "UOST" or len(set(ids)) != len(ids):
            raise ValueError("Event IDs must be unique strings, one per event")
        if not all(isinstance(value, str) and value for value in ids):
            raise ValueError("Event IDs must be nonempty strings")
        values = {}
        for name, descriptor in event_fields.items():
            is_grid = np.asarray(descriptor.values).ndim == 4
            dims = ("event", *SPATIAL_DIMS) if is_grid else ("event",)
            expected = (len(ids), *shape) if is_grid else (len(ids),)
            values.update(_data_vars({name: descriptor}, dims, expected))
            _physical_domain(name, values[name][1], active)
        ec = {"event": ("event", ids),
              "event_time": ("event", event_times, {"timezone": "UTC", "long_name": "exact discontinuity valid time"})}
        if any(len(value[0]) > 1 for value in values.values()):
            ec.update({axis: coords[axis] for axis in "xyz"})
        events = xr.Dataset(values, coords=ec, attrs={**attrs,
            "hydrostate:role": "ordered discontinuities; left/right values are never averaged",
            "start_datetime": _iso(event_times[0]), "end_datetime": _iso(event_times[-1])})
        if "x" in events.dims:
            events = spatial(events)
        nodes[EVENT_GROUP] = events
    if checkpoint_payload is not None:
        checkpoint = _jsonable(dict(checkpoint_payload))
        valid = _utc([checkpoint.get("valid_time_utc")], "checkpoint valid_time_utc")[0]
        if valid != times[-1]:
            raise ValueError("Checkpoint valid time must match the last model state")
        payload = json.dumps(checkpoint, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                             allow_nan=False)
        nodes[CHECKPOINT_GROUP] = xr.Dataset({"payload": ("checkpoint", np.asarray([payload]),
            {"units": "1", "media_type": "application/json", "long_name": "complete lossless model checkpoint"})},
            attrs={**attrs, "hydrostate:role": "complete restart; summed damage is not sufficient memory",
                   "hydrostate:checkpoint-protocol": checkpoint.get("schema", checkpoint.get("schema_version", "caller-defined")),
                   "start_datetime": _iso(valid), "end_datetime": _iso(valid)})
    tree = xr.DataTree.from_dict(nodes)
    validate_tree(tree)
    return tree


def validate_tree(tree: xr.DataTree) -> None:
    """Reject reference, grid, time, epoch and typed-variable incompatibilities."""
    if tree.attrs.get("hydrostate:contract") != CONTRACT or "/" + STATE_GROUP not in tree.groups:
        raise ValueError("Unsupported model datahub contract or missing model states")
    state = tree[STATE_GROUP].to_dataset()
    ref, crs = _reference(state.attrs.get("hydrostate:reference", {}))
    expected_grid_hash = _hash([state.x.values, state.y.values, state.z.values, state.active.values,
        {axis: float(state.attrs[f"hydrostate:d{axis}_m"]) for axis in "xyz"}])
    if state.attrs.get("hydrostate:grid-fingerprint") != expected_grid_hash:
        raise ValueError("State grid content disagrees with its fingerprint")
    for node in tree.subtree:
        ds = node.to_dataset()
        if not ds.sizes:
            continue
        _source_attrs(ds.attrs)
        if ds.attrs.get("hydrostate:grid-fingerprint") != state.attrs.get("hydrostate:grid-fingerprint"):
            raise ValueError(f"Grid fingerprint mismatch in {node.path}")
        if ds.attrs.get("hydrostate:reference-fingerprint") != _hash([ref]):
            raise ValueError(f"Reference fingerprint mismatch in {node.path}")
        if ds.attrs.get("hydrostate:reference") != ref:
            raise ValueError(f"Coordinate reference mismatch in {node.path}")
        if "time" in ds.dims:
            t = ds.time.values
            if t.dtype.kind != "M" or np.isnat(t).any() or ds.time.attrs.get("timezone") != "UTC":
                raise ValueError(f"{node.path} requires CF datetime64 with UTC metadata")
            if len(t) > 1 and np.any(np.diff(t) <= np.timedelta64(0, "ns")):
                raise ValueError(f"{node.path} time must increase uniquely")
        if "event" in ds.dims:
            ids = ds.event.values
            if len(set(ids)) != len(ids) or not all(isinstance(value, str) and value for value in ids):
                raise ValueError("Events require unique text IDs")
            et = ds.event_time.values
            if et.dtype.kind != "M" or np.isnat(et).any() or ds.event_time.attrs.get("timezone") != "UTC":
                raise ValueError("Events require exact datetime64 valid times with UTC metadata")
            if len(et) > 1 and np.any(np.diff(et) < np.timedelta64(0, "ns")):
                raise ValueError("Events must be ordered")
        if set(SPATIAL_DIMS) <= set(ds.dims):
            for axis in SPATIAL_DIMS:
                if not np.array_equal(ds[axis], state[axis]):
                    raise ValueError(f"{node.path} {axis} grid differs; remap explicitly upstream")
                spacing = ds.attrs.get(f"hydrostate:d{axis}_m")
                if (type(spacing) not in (int, float) or not np.isfinite(spacing) or spacing <= 0
                        or spacing != state.attrs.get(f"hydrostate:d{axis}_m")):
                    raise ValueError(f"{node.path} {axis} spacing must match the shared physical grid")
            if ds.active.dims != SPATIAL_DIMS or ds.active.dtype.kind != "b" or not np.array_equal(ds.active, state.active):
                raise ValueError(f"{node.path} active mask must match and remain boolean")
            if crs is None and ("spatial_ref" in ds or "proj:wkt2" in ds.attrs):
                raise ValueError("Synthetic local coordinates must not claim an earth CRS")
            if crs is not None:
                import rioxarray  # noqa: F401
                if ds.rio.crs is None or ds.rio.crs != crs:
                    raise ValueError(f"{node.path} requires matching spatial_ref CRS")
        for name, var in ds.data_vars.items():
            if name == "time_bounds":
                if var.dims != ("time", "nv") or var.shape[1] != 2 or var.dtype.kind != "M":
                    raise ValueError("time_bounds requires CF datetime64 (time,nv=2)")
                if np.isnat(var).any() or np.any(var.values[:, 1] <= var.values[:, 0]):
                    raise ValueError("Invalid interval time bounds")
                if (len(var) > 1 and np.any(var.values[1:, 0] < var.values[:-1, 1])):
                    raise ValueError("Interval time bounds overlap")
                if ds.time.attrs.get("bounds") != "time_bounds" or not np.array_equal(ds.time.values, var.values[:, 1]):
                    raise ValueError("Interval labels must identify their explicit ending times")
                continue
            if name == "spatial_ref":
                continue
            if name in FIELD_UNITS and var.attrs.get("units") != FIELD_UNITS[name]:
                raise ValueError(f"{name} has incompatible units")
            if set(SPATIAL_DIMS) <= set(var.dims):
                if var.dims[-3:] != SPATIAL_DIMS:
                    raise ValueError(f"{name} spatial dimension order must be z,y,x")
                _physical_domain(name, var.values, state.active.values)
            if var.dtype.kind not in "fiuUbOST":
                raise ValueError(f"{name} dtype cannot be coupled")
            if var.dtype.kind == "O" and not all(isinstance(value, str) for value in var.values.flat):
                raise ValueError(f"{name} contains nontext object values")
            if node.path == "/" + INTERVAL_GROUP and name != "active" and var.attrs.get("cell_methods") not in ("time: mean", "time: sum"):
                raise ValueError(f"{name} requires interval mean/sum metadata")
    if "/" + MATERIAL_GROUP in tree.groups:
        material = tree[MATERIAL_GROUP].to_dataset()
        for name in ("time", "material_epoch", "material_fingerprint"):
            if not np.array_equal(state[name], material[name]):
                raise ValueError(f"State/material {name} mismatch")
    static = tree[STATIC_GROUP].to_dataset() if "/" + STATIC_GROUP in tree.groups else xr.Dataset()
    dynamic = tree[MATERIAL_GROUP].to_dataset() if "/" + MATERIAL_GROUP in tree.groups else xr.Dataset()
    for i in range(len(state.time)):
        parts = []
        for ds, timed in ((static, False), (dynamic, True)):
            for name in sorted(ds.data_vars):
                if name in ("active", "spatial_ref"):
                    continue
                var = ds[name]
                parts.extend([name, var.attrs["units"], var.values[i] if timed else var.values])
        if state.material_fingerprint.values[i] != _hash(parts):
            raise ValueError("Material data disagree with the state's material fingerprint")
    if "/" + CHECKPOINT_GROUP in tree.groups:
        checkpoint = checkpoint_from_tree(tree)
        valid = _utc([checkpoint.get("valid_time_utc")], "checkpoint valid_time_utc")[0]
        if valid != state.time.values[-1]:
            raise ValueError("Checkpoint valid time must match the last model state")


def checkpoint_from_tree(tree: xr.DataTree) -> dict:
    """Return full caller-owned restart memory; caller validates its protocol/hash."""
    if "/" + CHECKPOINT_GROUP not in tree.groups:
        raise ValueError("The tree contains no complete model checkpoint")
    ds = tree[CHECKPOINT_GROUP].to_dataset()
    if ds.payload.dims != ("checkpoint",) or ds.payload.size != 1:
        raise ValueError("Checkpoint payload must contain exactly one JSON document")
    payload = ds.payload.values[0]
    if not isinstance(payload, str):
        raise ValueError("Checkpoint payload must be UTF-8 text")
    try:
        result = json.loads(payload)
    except (ValueError, TypeError) as exc:
        raise ValueError("Invalid checkpoint JSON") from exc
    if not isinstance(result, dict):
        raise ValueError("Checkpoint JSON must be an object")
    return _jsonable(result)


def save_tree(tree: xr.DataTree, path: str, *, overwrite: bool = False) -> None:
    """Validate, then use actual Gaia v3/consolidated/STAC storage APIs."""
    validate_tree(tree)
    try:
        from gaia_cli import io
    except ImportError as exc:
        raise ImportError("Install the pinned gaia-cli interoperability extra to write Gaia stores") from exc
    io.check_output(str(path), overwrite=overwrite)
    io.save_datatree(tree, str(path), overwrite=overwrite)


def open_tree(path: str) -> xr.DataTree:
    """Read with Gaia, then check model-level constraints Gaia does not enforce."""
    try:
        from gaia_cli import io
    except ImportError as exc:
        raise ImportError("Install the pinned gaia-cli interoperability extra to read Gaia stores") from exc
    tree = io.open_datatree(str(path))
    validate_tree(tree)
    return tree
