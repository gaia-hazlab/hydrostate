# Fixed-grid material changes, conserved water and hydraulic restarts

Hydrostate now has an executable interchange in
[`hydrostate.coupling`](../../src/hydrostate/coupling.py). Mechanicalstate can
replace independently prescribed material parameters while Hydrostate retains
ownership of hydraulic water, pressure and flux. This interface supports
synthetic evolution experiments. It does not solve mechanical deformation,
damage initiation, cyclic excess pore pressure, heat, vegetation recovery or
large landslide motion. Those processes need additional equations, evidence
and independently owned state.

## What is conserved

The existing mixed Richards backend integrates the dimensionless storage

\[
W(\psi; a)=\theta(\psi;a)+S_s\max(\psi,0),\qquad
H=z+\psi,\qquad \mathbf q=-K(\psi;a)\nabla H.
\]

Here \(a\) denotes the material parameters. For one constant-volume cell,
\(V\Delta W+\Delta t\sum Q_{\mathrm{out}}-\Delta t Q_{\mathrm{source}}=0\).
Physical liquid content is \(\theta\), physical saturation is
\(S_w=\theta/\phi\), and the elastic correction in \(W\) belongs in the
hydraulic ledger. It must not be substituted for \(\theta\) in a phase-volume
density or saturation closure. See the [flow derivation](flow_mathematics.md)
for the existing retention, conductivity and finite-volume equations.

At a prescribed instantaneous material change, no liquid crosses a face and
no water is injected by the change itself. The fixed-volume reconciliation is

\[
W_i(\psi_i^+;a_i^+)=W_i(\psi_i^-;a_i^-)=W_i^-\quad\text{for every active cell }i.
\]

Keeping \(\psi^-\) while changing retention or storage parameters generally
changes \(W\), inventing or deleting water. Keeping \(\theta\) instead loses
the elastic term. An updated pressure is therefore determined by the **new**
constitutive law and the explicitly retained **old** stored water. A
conductivity-only change leaves the storage law unchanged and preserves the
pressure exactly, including a saturated incompressible pressure reference.

These equations hold on the unchanged model bulk volume. They do not supply
the solid-mass balance, volumetric strain, stress-work balance or geometry of
a compaction/slide calculation. A prescribed porosity or retention change in
this interface is an external material scenario; it is not a prediction of
the mechanics that caused it. Changed cell coordinates, volumes or active
masks are rejected. A future moving-grid/remeshing interface must conserve
integrated water on intersecting volumes and account for deforming pore
volume, solid mass, material advection and relevant boundaries.

## Storage inverse and failure branches

For finite unsaturated pressure, \(\theta_r<W<\theta_s\). Set
\(m=1-1/n\) and \(S_e=(W-\theta_r)/(\theta_s-\theta_r)\). The exact equilibrium
inverse is

\[
\psi=-\alpha^{-1}\left[S_e^{-1/m}-1\right]^{1/n}<0.
\]

The implementation evaluates the bracket with `expm1(-log(Se)/m)` to retain
precision near saturation. For \(W\ge\theta_s\) and \(S_s>0\),

\[
\psi=(W-\theta_s)/S_s\ge0.
\]

Three states must not be hidden with clipping:

| Proposed storage | Consequence |
| --- | --- |
| \(W\le\theta_r\) | No finite pressure represents the requested state under this equilibrium law. Reject. |
| \(S_s=0, W>\theta_s\) | Storage exceeds the model capacity. Reject. |
| \(S_s=0, W=\theta_s\) | Any nonnegative pressure gives the same stored water. Require a declared nonnegative pressure hint. |

`reconcile_material_change` uses the old pressure as that explicit default
hint. A negative hint cannot establish a new saturated pressure; the update
fails unless the caller supplies an admissible nonnegative reference. The
linear manufactured-solution mode has the separate inverse \(\psi=W/C\).
Physical water and linear verification storage modes cannot be mixed.

Shrinking \(\theta_s\) when \(S_s\) is small can require thousands of metres of
pressure despite a mathematically valid inverse. Every operation checks a
declared finite pressure envelope. The default `(-100000, 100)` m is a
synthetic numerical guard, **not** a calibrated geological validity range.
Applications must choose an evidence-based pressure envelope and reject a
proposal that would require missing mechanics, boundary discharge or air
dynamics. Raising the bound is not physical validation. No rejected state is
clipped or committed, and neither the old model nor its pressure is mutated.

## Executable material transaction

```python
from hydrostate.coupling import reconcile_material_change
from hydrostate.flow import RichardsModel

# Build a replacement; never mutate old_model.theta_s in place.
candidate = RichardsModel(old_model.grid, **new_material_arrays)
update = reconcile_material_change(
    old_model, psi, candidate,
    stored_water=explicit_old_W,
    pressure_bounds_m=(-1000.0, 20.0),
)
# Install only after the complete transaction passes.
model, psi = candidate, update.pressure_head_m
```

The candidate receives actual full-grid material values, preserving exact
coordinates, spacings and active cells. Reference/datum identity must also be
checked by the enclosing typed spatial interface because `RichardsModel`
itself has no CRS or datum attribute. Neither pressure interpolation nor a
matching shape establishes grid/reference compatibility.

The result reports arrays `storage_before`, `storage_after`, `delta_storage`
and `delta_theta`, plus the signed sum, sum of absolute values and maximum
absolute **per-cell water-volume** residual. These prevent cancellation
between cells from hiding a nonconservative update. The inverse also verifies
the actual backend storage round trip to `5e-13` in dimensionless storage
per cell. This is a floating-point reconciliation tolerance, separate from a
flow solver's local volume-rate and pressure convergence tolerances. Pressure
and material hashes can be used by the enclosing state owner to prevent a
stale or mutated model from advancing.

## Hydraulic checkpoint, identity and clock

`build_checkpoint` returns detached standard JSON with protocol
`hydrostate.hydraulic-checkpoint/1`. It records:

- Exact grid, active mask, array order `zyx`, reference and SI/up conventions.
- Every actual backend material array, including scalar values expanded onto
  the domain and the optional mathematical linear-capacity branch.
- Pressure, physical water content and conserved \(W\), each with explicit
  units; inactive values are JSON `null`.
- UTC valid and issue times, the declared pressure envelope, metadata and
  grid/model/content SHA-256 fingerprints.

An issue time can precede a forecast valid time; these are different clocks.
Naive timestamps are rejected and aware offsets are normalized to UTC. The
hydraulic integration clock supports microsecond resolution; ISO fractions
beyond six digits and datetime subclasses with nonzero nanoseconds are rejected
instead of silently truncated. Generic DataHub files can retain nanoseconds,
but that storage precision does not change the model integration clock. A
state is not interpolated in time during restart. `restore_checkpoint` rebuilds
an actual backend, verifies units, geometry, material admissibility, content
integrity, and consistency of all three state arrays with the retention/storage
law. Optional `expected_reference`, `expected_grid`, `expected_model` and
`expected_valid_time_utc` checks reject a known target mismatch. Pass those
checks when restoring into an existing domain. The content hash detects
accidental changes, not authenticity or authorization.

```python
from hydrostate.coupling import build_checkpoint, restore_checkpoint

checkpoint = build_checkpoint(
    model, psi, valid_time_utc="2026-01-01T12:00:00Z",
    reference={"horizontal_crs": "LOCAL:synthetic", "vertical_datum": "local test datum",
               "length_unit": "m", "z_positive": "up"},
    metadata={"forcing_version": "synthetic-v1"},
)
restored = restore_checkpoint(checkpoint, expected_model=model)
```

This is a complete restart for the current equilibrium, history-free Richards
component. It is **not** a complete DHSVM–MODFLOW–mechanics restart. The
enclosing coupled restart must additionally save event IDs/order, damage and
healing memory, forcing position/version, exchange/analysis water ledgers,
boundary schedules, and every other process owner's physical state. Events
must be applied at their exact valid times between flow subintervals, and
restored events must not be applied twice. A seconds-long earthquake jump is
an external event here; the flow integrator does not resolve earthquake
shaking or cyclic pore-pressure generation.

## Verification performed

Run:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_coupling.py -v
```

The 15 focused tests pass using the actual backend. They check analytical
storage inverses across dry/unsaturated/saturated states, the incompressible
pressure ambiguity, unrepresentable storage and pressure-envelope rejection,
exact conductivity-only preservation, local water conservation under changed
retention and elastic storage, atomic rollback, material-mutation detection,
geometry rejection, and repeated material evolution without accumulated water.
Checkpoint tests exercise units, datums, clocks, masks, detached arrays and
constitutive consistency. An infiltration run with a changed retention law
produces exactly identical final pressure arrays and flux/storage ledgers
after checkpoint/restart compared with uninterrupted advancement using the
same substeps and boundary forcing.

These are software and constitutive checks. They establish neither the
physical correctness of a damage/healing closure nor transient nonlinear
solution mesh convergence, field calibration, hazard triggering or operational
forecast skill. Coupling-interval sensitivity, event ordering and damage-memory
restart belong in the enclosing temporal experiment; the existing numerical
and field-validation limitations continue to apply.
