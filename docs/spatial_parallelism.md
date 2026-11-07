# Spatial parallelism for a coupled vadose–groundwater model

Design assessment: 7 November 2026. Initial target: one workstation. The
schemes below are proposed implementation and benchmark choices, not an
implemented DHSVM–MODFLOW coupler or measured speedup. Spatial parallelism
divides one physical realization among workers; ensemble parallelism runs
different realizations. Benchmark them separately and together.

## Match the parallel scheme to the equations

| Physical model | Spatial work distribution | Information that must cross workers |
| --- | --- | --- |
| Independent 1-D Richards columns | Batch complete vertical columns across processes, threads or SIMD lanes | Surface forcing and groundwater boundary data; no lateral vadose exchange is represented |
| 2-D Richards hillslope section | Divide the section horizontally, keeping its vertical profiles together initially | Neighbor pressure/state data and conservative cross-partition face fluxes |
| 3-D Richards soil/vadose domain | Horizontal MPI tiles containing whole vertical columns; threads within tiles if beneficial | Lateral halos, residual/Jacobian contributions, global solver reductions and coarse-grid correction |
| Multilayer saturated aquifer | Connected MODFLOW models distributed using supported GWF–GWF exchanges and the extended parallel build | Head-dependent exchange terms in the coupled groundwater solution |

Independent columns are a **physical approximation**, not a general spatial
parallelization of 3-D Richards flow. They omit lateral unsaturated flow and
lateral redistribution of perched water. Use them only where observations and
comparison with a lateral-flow reference support that approximation. They may
be useful for vertical root-zone moisture, uptake and recharge, but that does
not establish slope-pressure or perched-water skill.

A **column-block preconditioner** is different: it can approximately solve the
vertical part of a linear system while the complete residual/Jacobian still
contains lateral flow. With a converged global solve, this retains the selected
2-D/3-D equations. Do not describe a preconditioner as removal of lateral
physics, or independent column runs as a full multidimensional solution.

ParFlow provides a relevant reference implementation of Richards flow and
terrain-following formulations. Its documented Richards computing topology
retains the vertical direction on each rank (`Process.Topology.R = 1`). This
supports testing horizontal decomposition first; it is not a claim that every
Richards solver must use that topology.
[ParFlow equations](https://parflow.readthedocs.io/en/latest/models.html),
[ParFlow computing topology](https://parflow.readthedocs.io/en/latest/keys.html#computing-topology).

## Baseline decomposition and data ownership

Start with a Cartesian `P_x × P_y × 1` MPI partition of the subsurface mesh.
Each rank owns a disjoint set of horizontal cells, their full vertical stacks,
local state/material arrays, and a halo around the tile. Co-locate the
corresponding land-surface work where practical. Keep a separate ownership map
for channel reaches and groundwater cells; their grids need not match.

For illustration, a `128 × 128 × 32` mesh split `2 × 2 × 1` gives each rank
`64 × 64 × 32` owned cells. A rectangular one-cell horizontal halo adds about
6.35% to the storage of a single cell-centered array on each rank. This is a
geometric example, not a memory estimate for a complete solver. Actual memory
also includes constitutive data, old/new states, Krylov vectors, matrices,
preconditioners and multigrid levels. A face-only stencil need not communicate
all corner entries in that rectangular allocation.

Keep vertical stacks on a rank initially because it avoids splitting local
vertical block solves and aligns surface forcing, root uptake and lower-boundary
exchange. It also avoids unnecessary vertical communication when there are
relatively few layers. Revisit the choice if very deep meshes, memory limits
or measured imbalance justify another decomposition; do not infer efficiency
from the rank count alone.

PETSc DMDA supplies structured-grid ownership and halo communication, including
7-point star and 27-point box stencils. It manages distributed arrays; it does
not supply the project's Richards equations, constitutive laws or conservative
interface implementation. Select halo width and corner exchange from the
actual discretization, including tensor/nonorthogonal terms.
[PETSc DMDACreate3d](https://petsc.org/release/manualpages/DMDA/DMDACreate3d/).

## Preserve face fluxes across tiles

Use a finite-volume residual with a documented storage law. One convenient
integrated convention is

```text
R_i = M_i(new) - M_i(old)
      + dt * sum(outward face flow rates)
      - dt * net external source rate = 0.
```

Here `M_i` is water storage in m³, including the chosen compressibility/storage
formulation; face flow rates are m³/s. For constant-density flow, evaluate
Darcy flux from total head `h = z + psi` and state-dependent conductivity.
Compute each shared physical face's signed flow once, or use an equivalent
deterministic evaluation that gives equal-and-opposite contributions to its
two cells. Assign a unique face ID and ownership rule. Summing all residuals
must cancel internal face transfers, including MPI and mixed-grid interfaces.

On a simple orthogonal, scalar-conductivity grid, resistance-weighted face
conductivity is a candidate. General anisotropic tensors, skewed meshes and
terrain-following coordinates need an appropriate consistent flux formulation;
do not assume a two-point harmonic rule is sufficient. Keep real elevation
and metric terms when following terrain. Tilting a grid without changing the
discrete gravity/geometry operators changes the physics.

During each residual/Jacobian evaluation:

1. Begin halo exchange of the current iterate and any required material/state
   data; compute owned interior terms while communication progresses.
2. Finish exchange, then compute boundary terms with the current neighbor
   values and assemble the paired interface contributions.
3. Complete global linear/nonlinear convergence checks. Exchange updated
   iterates as required by the solver; a single exchange per output hour is
   insufficient for an implicit multidimensional solution.
4. Accept the step only after local/global budgets and convergence criteria
   pass. A failed attempt restores the complete beginning-of-step state and
   retries with adjusted stepping/solver controls.

Do not independently solve each tile with frozen neighboring heads for an hour
and label the result a globally converged Richards solution. An iterated
domain-decomposition method can be valid, but needs interface convergence and
time-step sensitivity tests. If asynchronous local stepping is introduced,
integrate shared face volumes over a common interval and reconcile them before
acceptance; an unexplained mismatch is not numerical storage.

## Nonlinear and linear solution

Use a robust implicit serial reference before parallelizing. Compare modified
Picard and damped Newton for the chosen storage/retention laws, then use a
distributed Krylov solve with a suitable preconditioner. Monitor both storage/
flux residuals and pressure/head changes. Dry cells, wetting fronts, strong
layer contrasts and moving water tables can alter the Jacobian and iteration
cost substantially; reuse matrix sparsity, not stale physical coefficients.
PETSc supplies Newton with line search as an existing building block, not an
automatic Richards implementation.
[PETSc SNESNEWTONLS](https://petsc.org/release/manualpages/SNES/SNESNEWTONLS/).

Candidate preconditioners to compare are:

- Vertical line/column blocks plus a lateral/global coarse correction for
  thin, vertically refined domains. Block Jacobi alone may leave slowly
  varying lateral error unresolved.
- Geometric multigrid/semicoarsening for suitable logically structured grids;
  choose coarsening and smoothing from the actual anisotropy, grid spacing and
  boundary conditions, not only the labels horizontal and vertical.
- Algebraic multigrid for irregular connectivity or heterogeneous aquifers,
  checking setup memory, operator complexity and robustness across dry/wet
  states. More aggressive coarsening is a tuning experiment, not guaranteed
  improvement.

HYPRE documents SMG's line/plane smoothing and PFMG's cheaper point smoothing,
and BoomerAMG offers configurable parallel algebraic multigrid. These are
candidate **linear preconditioners**; they are not themselves complete
nonlinear vadose solvers. Their suitability for a particular nonsymmetric
Richards Jacobian must be demonstrated.
[HYPRE SMG/PFMG](https://hypre.readthedocs.io/en/latest/solvers-smg-pfmg.html),
[HYPRE BoomerAMG](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html).

For a more tightly coupled vadose–aquifer solver, consider a block Jacobian:

```text
[ J_vv  J_vg ] [delta psi_v] = -R_v
[ J_gv  J_gg ] [delta h_g  ]   -R_g
```

An approximate Schur complement can combine a vadose block solve with a
groundwater block solve while retaining the coupling terms. Apply approximate
inverse actions; do not build a dense Schur matrix for a watershed. PETSc
`PCFIELDSPLIT` supplies block/Schur infrastructure. Obtaining consistent
residuals, cross-derivatives and access to MODFLOW internals would be **new
integration/source work**, not a capability established by its array API.
Begin with conservative partitioned coupling and measure whether its iteration
cost justifies this more intrusive route.
[PETSc PCFIELDSPLIT](https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/).

## MODFLOW and the external vadose interface

MODFLOW 6's extended build supports spatial domain decomposition using MPI and
PETSc. Use supported GWF–GWF connections between aquifer submodels, supplying
cell IDs, connection geometry and the appropriate conductance formulation.
Use XT3D where justified for grid/interface geometry; compare partitioned and
unsplit solutions. The official parallel guide documents package limitations,
including GNC on parallel GWF–GWF exchanges, so pin and test the selected
release. A groundwater partition must not become an artificial no-flow
watershed boundary.
[MF6 extended guide](https://github.com/MODFLOW-ORG/modflow6/wiki/Extended-MODFLOW-User-Guide),
[GWF–GWF exchange](https://modflow6.readthedocs.io/en/latest/_mf6io/exg-gwfgwf.html).

A GWF–GWF exchange connects groundwater models; it does not automatically
connect an external Richards solver. Define that external interface separately:

- Every owned volume belongs to one solver. A Richards region can include
  locally saturated/perched cells, but do not duplicate their storage in
  DHSVM or MODFLOW. State what happens if the water table crosses the chosen
  interface or enters the root zone.
- Define intersecting faces/areas, elevation datums, pressure-to-head
  conversion, positive direction and the common exchange interval. Transfer
  integrated volumes conservatively; map boundary heads with their physical
  support rather than treating interpolation as a flux balance.
- Provide bidirectional head/flux feedback, interface storage treatment and
  coupling convergence. An imposed recharge series alone cannot transmit
  pressure or exfiltration feedback. Put ET and stream exchange under single
  process ownership.
- At each coupling iteration restore the beginning-of-interval state and
  reevaluate the trial interval. Commit states and the exchange ledger once,
  after convergence; never repeatedly advance or debit the same water.

Use persistent workers and in-memory arrays before replacing solvers. The
MODFLOW BMI exposes stepping and pointer access; current source documentation
warns that reinitialization in the same memory space is unsupported. Its XMI
documentation also says the `prepare_time_step(dt)` argument is ignored.
Verify those semantics in the pinned release and align actual TDIS/ATS times.
Test whether the desired parallel-library/API combination is supported rather
than assuming separate features compose automatically. DHSVM requires a
verified step/state interface and a global/static-state audit.
[MF6 BMI](https://modflow-org.github.io/modflow6/namespacemf6bmi.html),
[MF6 XMI](https://modflow-org.github.io/modflow6/namespacemf6xmi.html).

## Select physics spatially, with explicit interfaces

For a workstation, compare a regional groundwater backbone plus selected
Richards hillslope/valley patches with uniform 3-D Richards resolution. A 2-D
section is appropriate only where neglected cross-section flux is shown to be
unimportant. Use 1-D columns in areas where lateral vadose redistribution is
adequately small; use richer physics near layered perched systems, riparian
zones and geotechnical targets when observations support their parameterization.

Mixed physics requires care at patch edges. A saturated aquifer model cannot
alone supply an unsaturated lateral pressure boundary. A 3-D patch needs an
appropriate surrounding vadose solution, a sufficiently distant justified
boundary, or a validated interface approximation. If lateral flux is exchanged
between nominal columns, they are no longer independent column simulations:
the exchange must appear in the coupled residual and budget. One-way nests
do not establish two-way regional feedback.

Keep initial patch assignments fixed while validating the first implementation.
Changing physics or resolution during a run requires conservative transfer of
all relevant storage, consistent constitutive states and demonstrated absence
of artifacts. Do not switch a region to cheap physics simply because a hazard
threshold is about to be crossed.

## Load balance, memory and workstation execution

Partition by measured work as well as active-cell count. Snow/land-cover
options, root layers, wetting fronts, inactive terrain, locally refined layers,
streams and difficult groundwater connections can create unequal costs. Begin
with compact tiles and record per-rank work/wait times. If needed, weight
partitions using observed iteration/kernel costs, while limiting cut faces.
Repartition only after accounting for state migration, new maps and solver
setup costs. A moving wet front does not automatically justify dynamic
repartitioning every step.

Initially choose a common accepted timestep using the restrictions from the
coupled components, ending exactly at forcing/assimilation/exchange times.
The most demanding tile can limit a synchronized run. Local adaptive steps or
slower deep-aquifer updates are later experiments requiring interface-volume
consistency and error tests; they do not permit silently freezing groundwater.

Use contiguous state/work arrays and persistent buffers. For independent
columns, test batched tridiagonal kernels and layouts supporting both column
solves and vectorization across columns. For full multidimensional assembly,
test structure-of-arrays or blocked layouts against the actual stencil and
preconditioner. Avoid per-cell heap allocations and repeated Python loops.
Share only immutable forcing/geometry, and retain full restart/rollback state.
Budget multigrid/Krylov workspaces, not just the pressure array; matrix-free
Jacobian products can save matrix storage but still need an effective
preconditioner and verified residual evaluations.

Compare MPI-only, threads-only where safely supported, and a few MPI ranks with
threads. No worker count is prescribed without the actual workstation's RAM,
cores, memory bandwidth and library behavior. Within one realization, prefer
thread-private residual contributions or owned-face/gather assembly; shared
neighbor updates must not race. Across realizations, use separate communicators
and isolated mutable model state. Keep
`ensemble workers × spatial ranks per member × threads per rank` within the
tested resource allocation, including numerical-library threads.

GPU batching is a candidate for sufficiently large, regular workloads. Full
Richards execution also needs constitutive evaluation, residual/Jacobian work,
halos, solver/preconditioner support and acceptable transfer/memory costs.
PETSc GPU backends do not establish a GPU-ready DHSVM/MODFLOW coupler. Verify
the actual hardware/backend and pinned library compatibility; unified memory
does not remove capacity or bandwidth constraints.
[PETSc GPU support](https://petsc.org/release/overview/gpu_roadmap/).

## Benchmark and acceptance plan

Use identical equations, forcing, mesh, tolerances and outputs when comparing
parallel layouts. Changes to mesh/physics answer a separate accuracy–cost
question. Include hydrostatic equilibrium, vertical infiltration, an oblique
wetting front/layered hillslope, perched saturation, a moving water table and
gaining/losing stream exchange. Put material contrasts and active flow across
tile boundaries, not only within tiles.

| Experiment | Hold fixed | Vary | Report |
| --- | --- | --- | --- |
| Strong scaling | One mesh, physical problem and time interval | 1, 2, 4 and higher feasible spatial worker counts/layouts | Wall time, `T1/Tp`, `(T1/Tp)/p`, peak/total RAM, assembly/solve/communication/setup split, iteration counts |
| Weak scaling | Resolution, layers and approximately equal owned work per worker | Domain size and worker count together | Time per accepted step, imbalance, coarse-grid/reduction costs, iteration growth and RAM per worker |
| Ensemble–spatial allocation | Ensemble size, each member's physics and total resource budget | Many one-rank members versus fewer multi-rank members; thread settings | Completed members per hour, end-to-end cycle latency, RAM and forecast/DA cost |
| Selective physics | Defined prediction/validation targets | Columns, 2-D sections, 3-D patches and uniform reference | Pressure/moisture/recharge error, interface error, hazard-relevant state differences and cost |

For weak scaling, define a reproducible family of problems: enlarging a
physical domain can change conditioning, so equal cells per rank alone does
not isolate communication efficiency. Report that limitation. If a serial
case does not fit RAM, use a stated feasible baseline for relative scaling and
do not label it one-core speedup.

Measure warm-start and recurring work separately from initialization. Report
p50/p95/p99 end-to-end hourly-cycle times, failed attempts/recovery, data
ingestion, DA, output and delayed-data replay in addition to kernel timings.
Require local and global water closure, equal-and-opposite interface volumes,
restart equivalence and defined numerical agreement with the reference.
Scaling that increases failures or changes the accepted physics is not a
performance success. No runtime or memory reduction is established until
these measurements exist.
