# Vadose, ecological and geotechnical coupling design

Design date: 7 November 2026. This is the target architecture and its validation
plan. It is not evidence that the repository already implements these features,
that a toy run is calibrated, or that hourly operation has been demonstrated.
The companion [feasibility assessment](coupling_feasibility.md) discusses domain
sizes, observations and operational constraints.

## Physical target and implementation tiers

Use multidimensional, mixed-form Richards flow as the physical reference for
soil, vadose and shallow saturated dynamics. Three-dimensional flow is a
first-class requirement when lateral gradients, layered/perched saturation,
stream interactions or hazard-relevant pressure distributions require it. A
two-dimensional cross-section is a useful demonstrator and benchmark, with an
explicit out-of-plane width and boundary assumptions. Independent vertical
columns are a selectable approximation for domains where lateral vadose flow
has been shown to be unimportant; they are not the universal architecture.

The tiers below describe intended implementation scope, not interchangeable
claims of accuracy. Select processes from evidence and the quantity being
predicted, and record disabled processes with their justification.

| Tier | Processes and intended use | Acceptance before promotion |
| --- | --- | --- |
| A: transparent demonstrator | Isothermal, single-liquid-continuum 2-D Richards; heterogeneous retention/conductivity; rainfall and prescribed melt; ponding; root uptake with dry/wet stress; conservative face fluxes; pressure/suction and storage diagnostics | Water closure, hydrostatic equilibrium, mesh/time refinement and independent reference comparisons. Prescribed melt is not a snow energy model. |
| B: coupled research system | 3-D Richards with appropriate anisotropy/terrain geometry; DHSVM canopy/snow/energy interface; deep MODFLOW interface; stream/surface routing ownership; ecological compensation; complete restart and DA state | Interface and coupled budgets, converged exchange, water-table crossings, gaining/losing reaches, profile/head/streamflow validation and measured end-to-end throughput |
| C: evidence-triggered extensions | Root hydraulic redistribution/plant storage; explicit oxygen/solute or heat/vapor/ice; hysteresis; macropore/fracture domains; deformation feedback | Process-specific observations and benchmarks; parameter identifiability; added mass/energy balances and restart state |
| D: event and hazard modules | Stability, earthquake triggering/deformation, flood hydraulics, debris-flow initiation/runout or fire behavior, selected for the actual claim | Independent validation of the hazard quantity, with local geometry/materials and event-scale resolution |

Full 3-D variably saturated backends already exist. ParFlow is a primary
comparison for parallel integrated Richards/surface flow and terrain handling;
its large-cluster results do not establish workstation runtime for this system.
[Maxwell, 2013](https://doi.org/10.1016/j.advwatres.2012.10.001).
An October 2026 implementation also adds full Richards flow within MODFLOW 6,
including capillarity and lateral redistribution. The authors explicitly
identify it as research code outside official USGS releases. Evaluate it as a
candidate backend/reference, pin its source, and test compatibility with the
required API and packages; do not imply official support.
[Vázquez-Gasty et al., 2026](https://doi.org/10.1111/gwat.70112),
[implementation and benchmarks](https://www.gsienv.com/software/modflow-6-richards/modflow6-richards/).

## Conservation law, constitutive closure and discretization

Use SI units, elevation `z` positive upward, pressure head `psi` in meters and
total hydraulic head `H = z + psi`. For the rigid, constant-density,
isothermal liquid-water baseline with atmospheric pore air:

```text
q = -K_s(x) k_r(psi) grad(H)                  [m/s; Darcy flux vector]
partial theta(psi)/partial t + div(q) = f - S_r
theta                                        [m3 water / m3 bulk soil]
f, S_r                                       [1/s; external source, root sink]
u_w - u_a = rho_w g psi                       [Pa; u_a = 0 gauge in baseline]
```

`K_s` may be a tensor. Use a dimensionally consistent conductivity/retention
parameterization; saturated hydraulic conductivity is not intrinsic
permeability. If compressibility is retained, add an explicitly defined
storage contribution to the mixed balance and its water ledger. Do not add
MODFLOW specific yield to the same upper-domain water-table drainage already
represented by the retention curve. A saturated incompressible region has an
elliptic pressure equation; forcing a positive arbitrary moisture capacity to
avoid that limit changes the physical storage.

For an initial single-domain van Genuchten–Mualem closure:

```text
S_e = (theta - theta_r) / (theta_s - theta_r)
S_e(psi < 0) = [1 + (alpha |psi|)^n]^(-m), m = 1 - 1/n
S_e(psi >= 0) = 1
theta = theta_r + (theta_s - theta_r) S_e
k_r = S_e^ell [1 - (1 - S_e^(1/m))^m]^2
alpha > 0 [1/m]; n > 1; 0 <= theta_r < theta_s <= porosity
```

`ell = 0.5` is a conventional initial choice, not a measured universal value.
Effective saturation `S_e` is not the liquid pore saturation `theta / porosity`;
export both with their definitions when ecology or mechanics consumes them.
Fit/justify retention and conductivity using material data and uncertainty;
texture-class lookups alone are weak evidence for volcanic layers, ash,
fractures or engineered fill. Stabilize evaluations near saturation/dryness
without silently clipping physical states or adding substantial artificial
conductivity. [van Genuchten, 1980](https://doi.org/10.2136/sssaj1980.03615995004400050002x),
[Mualem, 1976](https://doi.org/10.1029/WR012i003p00513).

Discretize cell-integrated water storage and one shared flux per internal
face. For cell volume `V_i` and outward face flux `Q_if`:

```text
V_i [theta_i(new) - theta_i(old)]
  + dt sum_faces Q_if(new) - dt V_i [f_i(new) - S_ri(new)] = 0
```

This mixed residual uses the actual change of water content, rather than
replacing it by `C(psi_old) delta_psi`. Use backward Euler initially with
modified Picard or damped Newton iterations, adaptive step rejection/retry and
separate storage/flux and pressure convergence checks. Mixed form supports
conservation, but finite solver tolerances and boundary treatments still need
measured residuals. [Celia et al., 1990](https://doi.org/10.1029/WR026i007p01483).

Compute flux from **total head**, preserving hydrostatic no-flow equilibrium.
Document face conductivity averaging and test sharp contrasts and dry fronts.
Orthogonal, aligned scalar/diagonal-conductivity toy grids can use a two-point
flux stencil. Nonorthogonal terrain grids or rotated anisotropic tensors need
a consistent discretization; a simple two-point stencil is not generally
adequate. Soil depth coordinates on a slope must retain the gravity projection
and metric terms. Refine vertically at horizons, roots, slip surfaces and
water-table transitions, rather than assigning equal thickness everywhere.

## Surface, snow and ecological processes

**Canopy and snow.** DHSVM owns interception, precipitation phase, snow mass
and energy, canopy evaporation/sublimation and melt production when coupled.
It sends throughfall, liquid melt and energy-consistent remaining soil
evaporation/transpiration demand, rather than raw precipitation plus a second
full PET allocation. Glacier melt/routing is included where material. A toy
prescribed melt pulse tests liquid delivery only. Full snow ownership remains
with one component. [PNNL DHSVM](https://github.com/pnnl/DHSVM-PNNL).

**Surface boundary.** Represent a nonnegative ponded depth and use a switching
atmospheric-flux/ponded-pressure condition. Dry-surface evaporation must reduce
when the near-surface soil cannot supply the demand. Rainfall exceeding
infiltration enters ponding and then explicit runoff; positive subsurface
pressure can cause exfiltration. Exchange the same signed top-face volume in
both subsurface and surface ledgers. For a surface patch, with downward
subsurface flux `f_top`, actual atmospheric evaporation `E_surface`, liquid
delivery `I` and runoff export `R`:

```text
dp/dt = I + runon - E_surface - f_top - R
```

When dry, the signed soil-top flux includes the actual soil evaporation;
`E_surface` records that same atmospheric removal so that the zero-storage
surface balance remains consistent. When ponded, evaporation is drawn from
the pond and `psi_surface = p`. Implement and test these boundary branches
together; do not remove soil evaporation again as an interior sink. Resolve
event forcing and boundary switches rather than smearing a short intense
burst into an hourly mean. A local overflow threshold in a toy model is not
channel/floodplain hydrodynamics.

**Root water uptake.** Keep potential transpiration distinct from bare-soil
evaporation and interception loss. For a plant group with volumetric potential
transpiration `Q_p [m3/s]`, use a normalized root distribution `b(x)` with
`integral b dV = 1`; a baseline sink is
`S_r = Q_p b(x) a_dry(psi) a_wet(psi)`. Define thresholds, units and their
dependence on demand/species. A Feddes-type wet-end response represents
anaerobiosis empirically; it is **not** an oxygen diffusion/respiration model.
Use soil/plant observations to identify dry/wet response and rooting depth.
[HYDRUS-1D technical manual, root uptake](https://pc-progress.com/Downloads/Pgm_Hydrus1D/HYDRUS1D-4.15.pdf).

**Compensation.** Allow wetter accessible roots to increase uptake within
root/soil hydraulic capacity, while total uptake remains bounded by demand
and available water. Normalizing stressed root weights to force full demand
even when all roots are stressed is not a physical compensation model. A
hydraulic-architecture model can represent uptake redistribution through
root-system and compensatory conductances; these require additional
identifiable parameters. [Couvreur et al., 2012](https://doi.org/10.5194/hess-16-2957-2012).

**Hydraulic redistribution.** A separate root/xylem module can calculate signed
soil–root exchange from water-potential differences. It may return water from
roots to dry soil; it is not a nonnegative transpiration sink. At zero
transpiration and negligible plant storage, redistribution sums to zero across
the connected roots. If plant storage is retained, the net soil removal equals
transpiration plus plant-storage increase. Retain nighttime transpiration when
supported. Root networks connecting patches also create spatial coupling and
must participate in the parallel exchange scheme.
[Amenu and Kumar, 2008](https://doi.org/10.5194/hess-12-55-2008).

**Ecological outputs and optional feedback.** Export actual/potential
transpiration, uptake by depth, root-zone water potential, stress duration,
waterlogging duration and available-water indicators. Prescribed seasonal
roots/LAI are a first implementation; dynamic phenology, growth, rooting,
stomatal/plant hydraulics, carbon/nutrient cycling and mortality require
additional ecological modules and calibration. Soil moisture is not measured
live/dead fuel moisture. Root cohesion is an independently parameterized
mechanical effect, not an automatic consequence of uptake.
Roots, biopores, litter and mortality can also change retention, conductivity,
anisotropy and mechanical reinforcement. Represent supported structural
feedback through evidenced, time-varying material/ecological parameters and
associated storage changes; changing an uptake sink alone does not capture it.

## Additional processes: trigger and contract

| Process | Evidence that should trigger inclusion | Added equations/state and constraints |
| --- | --- | --- |
| Heat and vapor | Soil temperature/evaporation or dry-surface moisture cannot be explained by liquid flow alone | Energy balance, thermal conductivity/capacity, liquid/vapor transport and surface energy fluxes; track latent heat and distinguish liquid from total water |
| Freeze/thaw | Seasonal frost or glacier-adjacent freezing controls infiltration, storage or hazard conditions | Liquid/ice partition, enthalpy/latent heat, freezing curve and ice-dependent conductivity; temperature alone is not an ice-storage model |
| Hysteresis/air entrapment | Wetting/drying cycles show different retention or infiltration/pressure response | Main/scanning curves, reversal history and entrapped-air state; restore history at restart and after DA |
| Macropores/preferential flow | Rapid bypass, tracer/profile evidence or fractures violate local single-continuum equilibrium | Dual-permeability or explicit pathways with exchange and separate stores; a single larger `K_s` is not a general substitute |
| Fracture/matrix exchange | Volcanic/bedrock structure controls nonlocal recharge, perched water or spring response | Evidence-based equivalent tensor, dual continuum or resolved fractures, with matrix exchange and appropriate boundaries |
| Oxygen/salinity | Root stress requires aeration or osmotic effects beyond pressure proxies | Gas/dissolved oxygen transport plus consumption, or solute transport and osmotic plant response; do not label wet-end stress a measured oxygen concentration |
| Density/viscosity changes and hydrothermal flow | Temperature/salinity/pressure materially changes buoyancy or fluid properties | Pressure-based mass flow with density/viscosity and energy/solute coupling; a fixed-density head equation is insufficient for buoyancy-driven circulation |
| Deformation and changing porosity | Compaction, shrink/swell, dilation or seismic loading changes storage/flow or stability | Consistent poromechanical mass/force balance and evolving geometry/retention/permeability; account for water in changing pore volume |
| Fire and recovery | Burn severity and measurements support changed hydrology/ecology | Time-varying canopy/litter, repellency, conductivity, roughness, erodibility and root reinforcement; do not assume all burns reduce conductivity |

These extensions are optional only when immaterial to the selected application.
SUTRA 4 provides a published reference for saturated/unsaturated
variable-density flow with energy and freeze/thaw; HYDRUS documents
nonequilibrium/hysteresis options. Reuse validated components where possible.
[USGS SUTRA 4 enhancements](https://doi.org/10.3133/tm6A63),
[HYDRUS-1D tutorial](https://www.pc-progress.com/Downloads/Public_Lib_H1D/HYDRUS-1D_Tutorial_V1.00_2018.pdf).
Richards flow assumes an adequately represented pore-air condition. If air
compression or trapped-gas dynamics controls infiltration, assess two-phase
flow rather than quietly retaining atmospheric air pressure.

## Nonoverlapping DHSVM–Richards–MODFLOW contract

The preferred multidimensional research partition is a **fixed spatial split**:
upper Richards owns soil, all vadose storage, perched water and any saturated
flow above the split; MODFLOW owns the disjoint deeper saturated volume. Choose
the split below the lowest expected regional water table, supported by
hydrogeology and spin-up. Water tables can cross upper Richards cells without
moving a grid interface or transferring their storage to MODFLOW.

At the common interface, enforce continuity of hydraulic head and equal,
oppositely signed normal water flux. Evaluate finite-volume connection fluxes
using the adjacent half-cell resistances and any evidenced contact resistance;
cell-center head equality alone is not the interface condition. Transfer
integrated volumes with conservative area mapping. Use a converged
Dirichlet–Neumann/Robin iteration, or a verified monolithic interface. Retain
beginning-of-step states and commit storage once after convergence.

This partition requires saturation at the split and in the lower GWF domain.
Monitor that condition. If the split becomes unsaturated, the assumption has
failed: deepen/redesign the split or adopt a full variably saturated backend
before accepting that regime. Do not silently represent the missing interval
by instantaneous recharge, or maintain overlapping Richards and MF6 storage.
A dynamically moving partition is a separate research implementation requiring
conservative state remapping, retention/specific-yield reconciliation and
crossing benchmarks; it is not the default.

| Component | Owned quantities | Interface and restrictions |
| --- | --- | --- |
| DHSVM land surface | Canopy, snow/ice mass/energy and assigned surface energy demand | Send liquid delivery and allocated evaporation/transpiration demand. Replace overlapping soil/saturated routing where Richards owns those volumes. |
| Richards domain | Root-zone/vadose/upper saturated water, pressure, internal lateral/vertical fluxes and assigned root sinks | Return actual ET, pressure/soil state and surface exchange; export deep interface volumes. Do not run UZF on the same storage. |
| Surface/channel owner | Ponding/overland/channel/lake storage and routing, as explicitly partitioned | Exchange runoff, exfiltration and streambed leakage once; MF6 SFR is an option if its routing assumptions suit the application. |
| MODFLOW GWF | Deep saturated heads, elastic storage, pumping and lower lateral/boundary flows | Receive/return conjugate interface volumes and heads; do not add specific yield for upper Richards water-table storage. |
| Ecological module | Assigned uptake/plant state, including storage if enabled | One demand allocation; signed redistribution exchanges cancel except transpiration/plant storage. |
| Mechanics/hazard modules | Strength/stress, deformation, damage, initiation or runout states | Hydrology supplies pressures/suction/saturation with support and uncertainty. Feedback is off until implemented and validated. |

Alternative shallow/deep partitions retaining DHSVM soil routing are allowed
only for physically distinct stores. If an intervening vadose module starts
below that soil, feed it soil-bottom drainage and allocated unmet uptake, not
raw rainfall or duplicate soil storage. Independent 1-D columns/UZF can be
benchmarked as faster choices for suitable regions. Standard UZF is
gravity-dominated vertical kinematic flow; it does not supply general
capillary-fringe or lateral unsaturated physics.
[USGS MF6 theory](https://pubs.usgs.gov/tm/06/a55/tm6a55.pdf).

The MF6 API supports controlled stepping and variable access, but allowed
writes, pointer lifetime, step control and multiple-instance behavior must be
verified for the pinned release. Proposed component methods are `initialize`,
`advance_trial`, `get_state`, `get_flux_integrals`, `commit`, `rollback`,
`checkpoint` and `finalize`; these names are our contract, not existing DHSVM
or MF6 API promises. Every exchange includes CRS/datum, geometry, units,
sign, valid-time window and ownership. [MODFLOW API](https://doi.org/10.1016/j.envsoft.2021.105257).

## Geotechnical products and optional feedback

Export coordinate/depth-resolved `psi`, `u_w`, saturation, water content,
porosity, hydraulic gradients and seepage fluxes, together with layer identity,
uncertainty and time support. Resolve candidate slip surfaces rather than
reconstructing all pressure from one regional water-table elevation. Distinguish
positive water pressure from suction and use a declared pore-air reference.

For saturated soils, an effective-stress calculation typically uses
`sigma' = sigma - u_w I` under the appropriate total-stress convention. For
unsaturated soils, use a justified stress/strength closure, for example a
calibrated Bishop-type weighting or suction-stress relation. Do not apply the
full negative water pressure as equivalent cohesion, or assume `chi = S_e`
universally. Retention data alone do not validate shear strength.
[Lu and Likos, 2006](https://doi.org/10.1061/(ASCE)1090-0241(2006)132:2(131)).

The initial contract is one-way hydrology to mechanics. A factor of safety is
a stability diagnostic, not displacement/runout. Liquefaction and lateral
spreading also need shaking, geotechnical resistance and event-scale mechanics;
Richards plus hourly groundwater does not generate cyclic excess pressure.
[USGS mechanisms](https://www.usgs.gov/programs/earthquake-hazards/what-are-effects-earthquakes),
[TRIGRS scope](https://pubs.usgs.gov/of/2008/1159/).
Two-way mechanics requires deformation-induced storage/pressure, porosity and
permeability changes with energy/momentum consistency and a tested coupling
scheme. Its event time steps are separate from the hourly antecedent cycle.

## State, observations, parallel execution and optimization

Checkpoint primary pressure and consistent water storage, ponding/surface
stores, snow/energy state, aquifer state, roots/plant storage and all activated
history variables. Include material parameters, mesh/version, forcing cursor,
accepted model time and pending exchange volumes. Hysteresis, ice, solutes and
root hydraulic state cannot be reconstructed from water-table head alone.

DA operators use sensor depth/volume, water-content measurement physics, well
screen/datum, routed discharge and forcing accumulation windows. Pressure and
moisture increments must satisfy constitutive laws; project constrained states
without hiding clipping losses. Record every analysis-induced water/energy
increment and preserve a replay window for delayed observations. Parameter
updates may change the retention relation and therefore require explicit state
reconciliation. Validate unobserved depth profiles and held-out sites/events,
not just fit at assimilated sensors.
[Mass-constrained assimilation](https://doi.org/10.5194/hess-16-105-2012).

For 3-D flow, partition spatially with ghost/halo exchange and one authoritative
face flux; reduce global residuals and synchronize accepted steps. Domain
decomposition needs an effective nonlinear/linear solver and preconditioner,
especially across large conductivity contrasts. Cache mesh topology/sparsity,
not state-dependent conductivity, storage or coupling conductance. Profile
assembly, solver iterations, interface iterations, memory traffic and I/O.
Independent columns can be batched/parallelized only in the declared column
approximation; turning off lateral flux is a change of equations, not a
parallelization technique. Root networks, runoff and channels can also connect
partitions.

Use isolated persistent ensemble workers where model global/static state is
not reentrant. Compare within-member spatial parallelism with member batching
on the actual workstation; cap numerical-library threads to avoid
oversubscription. Report serial/parallel numerical agreement, local/global
budgets, peak RAM and tail cycle times. Benchmark MF6 plus vadose plus
assimilation, rather than extrapolating from one kernel or one synthetic run.
See [spatial parallelism](spatial_parallelism.md) for detailed decomposition and
exchange schemes.

## Staged validation and decision gates

1. **Constitutive and balance checks:** invert retention consistently; recover
   saturated limits; verify hydrostatic equilibrium, equal/opposite internal
   fluxes, ponding/evaporation switching and cell/global budgets. Include zero
   demand and zero-input cases.
2. **Flow benchmarks:** compare infiltration/drainage and capillary rise with
   an independent Richards reference; test dry starts, fine-over-coarse and
   coarse-over-fine layers, low-permeability lenses, perched water and lateral
   redistribution. Repeat in 2-D/3-D as appropriate; a column benchmark alone
   does not validate multidimensional flow.
3. **Ecological checks:** validate drought/wet stress and uptake by depth;
   compensation cannot exceed demand/capacity. Test zero-demand root
   redistribution conservation and plant storage. Compare oxygen/thermal
   outputs only if those equations are actually enabled.
4. **Coupling checks:** use a nonoverlapping reference problem with a known
   interface solution; test rising/falling water tables, upward/downward
   exchange, split-saturation failure detection, strong storm feedback and
   stream gain/loss. Demonstrate mesh, time-step and coupling-interval
   convergence, rollback and restart equivalence.
5. **Field and DA checks:** spin up/calibrate offline; validate pressure,
   moisture profiles, heads, SWE and discharge across wet/dry seasons and
   extremes; withhold sites/events. Measure assimilation storage increments,
   posterior consistency and forecast reliability.
6. **Hazard and operational checks:** validate the particular triggering,
   deformation or inundation claim and its event resolution. Measure complete
   ingestion-to-product cycle cost on the target hardware, including replay,
   failures and forecasts. Promote a toy only as far as the evidence supports.

Set signed volume/depth and normalized budget tolerances for each experiment
before acceptance, with the normalization stated. Inspect maximum cell/interface
errors and cancellation as well as global closure. A faster solver is accepted
only with defined numerical/physical equivalence and demonstrated benefit.
