# Hydrostate–Mechanicalstate model data in Gaia DataHub

The implemented contract is `hydrostate.mechanicalstate.datahub/1` in
[`datahub.py`](../src/hydrostate/datahub.py). It uses Gaia's actual
`gaia_cli.io.save_datatree` / `open_datatree` APIs, rather than a similarly named
local writer. Compatibility was checked against gaia-cli commit
[`aa8ed3def5af4c32a4a1e601252fcee439344ac8`](https://github.com/gaia-hazlab/gaia-cli/tree/aa8ed3def5af4c32a4a1e601252fcee439344ac8).
This is a typed model-output and restart extension of Gaia storage conventions;
the groups are not registered `gaia stage` upstream datasets.

Gaia's source requirements are in
[adding a dataset](https://github.com/gaia-hazlab/gaia-cli/blob/aa8ed3def5af4c32a4a1e601252fcee439344ac8/dev/adding-a-dataset.md),
[provenance metadata](https://github.com/gaia-hazlab/gaia-cli/blob/aa8ed3def5af4c32a4a1e601252fcee439344ac8/dev/provenance-metadata.md)
and its actual
[I/O implementation](https://github.com/gaia-hazlab/gaia-cli/blob/aa8ed3def5af4c32a4a1e601252fcee439344ac8/src/gaia_cli/io.py).
The storage roadmap describes alternatives and older string behavior; executable
I/O code at the pinned commit takes precedence for this integration.

## DataTree layout and ownership

| Group | Dimensions | Meaning and owner |
| --- | --- | --- |
| `rasters/hydrostate-state` | `time,z,y,x` | Instantaneous pressure head, volumetric water content, total storage and hydraulic head; Hydrostate owns water evolution. |
| `rasters/mechanicalstate-materials` | `z,y,x` | Reference structural/hydraulic properties, facies and provenance. |
| `rasters/mechanicalstate-state` | `time,z,y,x` | Evolving hydraulic properties, damage and geophysical diagnostics at the same state valid times. |
| `timeseries/hydrostate-fluxes` | `time` or `time,z,y,x`; `time_bounds(time,nv)` | Interval-integrated water ledger or interval-mean rate, explicitly distinguished by units and `cell_methods`. |
| `events/material-jumps` | `event` or `event,z,y,x` | Exact discontinuity timestamps, unique IDs, epochs, signed water adjustment and separately named left/right states. |
| `checkpoints/hydro-mechanical` | one `checkpoint` record | A complete lossless JSON restart, including recovery modes and event history; not a reconstruction from displayed fields. |

Spatial leaves contain the same boolean `active(z,y,x)` mask. A single uniform,
ascending Cartesian grid is used by the present coupling. Cell centers, cell
widths and active geometry are fingerprinted. All groups carry the same grid
fingerprint, complete reference metadata and reference fingerprint. Material
fingerprints bind array values, units, field names and numeric types; the state
and evolving-material groups share per-time material epochs and fingerprints.
On reading, the adapter recomputes grid and material fingerprints and rejects
mismatches. A hash checks identity, not the accuracy or uniqueness of a
petrophysical interpretation.

The adapter does **no** coordinate transformation, vertical datum conversion,
spatial upscaling, state remapping or interpolation in time. A common DataTree
format alone is insufficient to justify any of these operations.

## Types, units and missing values

The interface requires explicit units on each field. Known quantities have a
strict unit registry; for example, pressure head is `m`, pore pressure `Pa`,
conductivity `m/s`, intrinsic permeability `m2`, velocity `m/s`, density
`kg/m3`, retention alpha `1/m` and volumetric water content `m3/m3`.
Physical total porosity is `1`. Saturation and porosity are different quantities.
The independent storage variable `total_water_storage_fraction` represents
`W = theta + Ss max(psi,0)` for the current kernel; it is not silently relabeled
as water content. Arbitrary field metadata cannot invent CF `standard_name`s;
use descriptive `long_name`, units, support and provenance. Time has the actual
CF standard name `time`.

Float32 observations/diagnostics and float64 solver fields retain their widths.
Registered prognostic and physical continuous material fields require finite
floating values in every active cell and NaN in every inactive cell; an active
pressure/storage hole is rejected instead of being passed to a solver. Unknown
analysis quantities may retain missing values with appropriate quality metadata.
NaN is the floating missing value; infinity is rejected. Booleans remain
booleans, including false inactive cells. Integer categories retain their integer
width and explicit sentinel values; their physical validity comes from the
`active` mask and optional `missing_value_code` metadata. Text is UTF-8 and
round-trips without fixed-width truncation. NumPy's native variable string dtype
(`kind='T'`), fixed-width Unicode and Python string object arrays are supported;
arbitrary object arrays and `None` values are rejected. This contract does not
zero-fill missing data.

**Integer fill semantics require care.** A canonical integer category is not a
CF-masked float: applying CF `_FillValue` decoding to an encoded integer is
allowed to return floating values with NaN. If raw encoded integer widths are
needed for a separately packed external product, read that product with
`mask_and_scale=False` and interpret its mask explicitly. At the pinned Gaia
commit, the small-array writer replaces same-dtype integer `fill_value` encoding
with its chunk encoding; a requested integer fill can therefore revert to the
Zarr codec default (zero). The canonical model category groups use a separate
active mask and explicit integer sentinel instead of depending on that fill.
Gaia writes complete empty/fill chunks, so codec fill is not the physical active
mask. The round-trip tests verify the actual decoded typed values. They do not
establish GDAL nodata behavior for three-dimensional model categories; readers
that treat codec zero as missing can misinterpret a legitimate category zero.

## Coordinates and time

Real georeferenced models require a parseable projected horizontal CRS with
metre axis units and an explicitly declared vertical datum. The adapter writes
CF `spatial_ref` with rioxarray. Gaia then writes horizontal GeoZarr `proj:` and
`spatial:` metadata and a STAC polygon/bounding box. Its horizontal GeoZarr shape
and transform describe `y,x`; they do not define a vertical transform.
`z` is elevation in metres, positive upward, with the vertical datum recorded.
Providing CRS metadata does not demonstrate that a datum conversion occurred.

Synthetic demonstrations use a reference such as
`horizontal_crs='LOCAL_SYNTHETIC_CARTESIAN'`, an explicit synthetic vertical
origin, `length_unit='m'` and `z_positive='up'`. These leaves have no `spatial_ref`
or earth CRS, and their STAC Item has no earth geometry. Declaring both a real
CRS and `synthetic_local=True` is rejected. Synthetic local coordinates must
never be assigned an EPSG code merely to make a map viewer accept them.

Timestamp inputs must be timezone-aware ISO text or datetime objects; naive
datetimes and bare datetime64 inputs are rejected. The adapter converts to UTC
and stores CF-decodable `datetime64[ns]` with `timezone='UTC'`. ISO nanoseconds are
preserved, including event timestamps between hourly outputs. Snapshot times
are unique and increasing. Coincident event times are permitted with unique
IDs and explicit event order.

Instantaneous states use `cell_methods='time: point'`. Interval fields require
explicit ordered nonoverlapping start/end bounds and `time: mean` for rates or
`time: sum` for integrated quantities. Their labels are the interval-ending
timestamps. Bounds may have gaps; a missing interval is not invented. Event
left/right states reside on the independent `event` axis and are not averaged
into hourly states. A complete checkpoint's `valid_time_utc` must match the
last state snapshot.

## Gaia alignment is not a coupling operation

**Do not run `gaia align` on these prognostic-state, material, event or restart
groups.** A `hydrostate:alignment-policy` attribute states this restriction;
current Gaia does not enforce that marker. It remains a caller obligation.

At the pinned commit, Gaia's temporal aggregation treats point samples as
representative interval-ending values, defaults to means, and casts nonfloating
timed variables to float32. Those operations can erase event jumps and epochs.
Gaia's regridding is horizontal. Its conservative option is a normalized
overlap-area mean that skips missing source values; it is not a conservative
three-dimensional transfer of cell water volumes, moving geometry or fluxes.
Some reprojection paths accept only two or three array dimensions. None supplies
vertical materials upscaling or a physically reconciled restart.

Use Gaia alignment for appropriate **derived forcing or comparison views**, with
their accumulation windows, coverage and category semantics checked. Coupling
must preserve `W`, interval volumes, face signs, process ownership and exact
event timing through model-aware operations. The kernel and coupling tests
establish those separate constraints for the implemented synthetic domain.

## Provenance, storage and complete restart

Every leaf has source `title`, `description`, `license`, `providers`, `links`
(license and `derived_from`/`via`), `processing:lineage` and time extent.
`source_attrs(...)` supplies a minimal synthetic-model template; actual input
licenses and sources must replace it for a field application. Model versions,
material provenance and assumptions should be included in source metadata and
the complete checkpoint. Gaia's writer adds actual `created` and
`gaia:cli-version` run metadata, writes Zarr v3 with consolidated metadata, and
writes the neighboring STAC Item. `hydrostate:gaia-contract-commit` records the
audited contract; `gaia:cli-version` records the writer actually used.

The checkpoint payload is one complete JSON object, stored as UTF-8 with media
type `application/json`. It preserves arbitrarily long text and all restart
memory declared by the caller. It must include individual healing/recovery
mode amplitudes and timescales: their sum does not determine future evolution.
`checkpoint_from_tree(tree)` returns that object; the owning coupled runner must
also validate its schema, content hash, material identity and physics invariants
before resuming. This adapter does not infer restart variables from pictures or
summary time series. Missing checkpoint data causes a restart error.

## Reproduction and tested scope

```python
from hydrostate.datahub import Field, build_tree, source_attrs, save_tree, open_tree

# grid is a validated Hydrostate Grid3D, and solver_arrays are [time,z,y,x].
tree = build_tree(
    grid=grid, reference=reference,
    provenance=source_attrs(
        title="Synthetic coupled valley", description="Research demonstration",
        source_commit=hydrostate_commit, lineage="Exact-grid hydromechanical experiment"),
    state_times_utc=valid_times,
    state_fields={"pressure_head_m": Field(psi, "m"),
                  "theta": Field(theta, "m3/m3"),
                  "total_water_storage_fraction": Field(W, "m3/m3")},
    static_fields=reference_materials,
    material_fields=evolving_materials,
    material_epochs=epochs,
    interval_bounds_utc=intervals,
    interval_fields={"boundary_out_volume_m3":
        Field(outflow_volumes, "m3", {"cell_methods": "time: sum"})},
    event_times_utc=event_times, event_ids=event_ids, event_fields=event_fields,
    checkpoint_payload=complete_restart)
save_tree(tree, "coupled.zarr")
restored_tree = open_tree("coupled.zarr")
```

Install the repository's `datahub` extra for xarray/Gaia runtime support. The
mandatory interoperability job uses the pinned Gaia source; setting
`HYDROSTATE_REQUIRE_GAIA=1` makes a missing Gaia import a test failure rather than
a skipped integration test. To use an already installed Gaia environment:

```bash
HYDROSTATE_REQUIRE_GAIA=1 PYTHONPATH=src \
  python -m unittest discover -s tests -p test_datahub.py -v
```

The tests write and reopen real local and EPSG:32610 stores through the actual
Gaia APIs. They check float widths, integer categories, false booleans, NaNs,
UTF-8 labels, exact event nanoseconds, full long Unicode restart payloads,
time bounds, shared material epochs/fingerprints, consolidated Zarr v3,
GeoZarr/CF metadata, provenance and STAC earth geometry only for real CRS.
Rejected cases include unit mistakes, naive timestamps, implicit feet/degree
coordinates, mismatched grids, mutated material arrays, invalid epochs,
overlapping intervals and stale checkpoints. These are format/identity tests,
not field validation, imaging inversion validation or a proof of hazard skill.
