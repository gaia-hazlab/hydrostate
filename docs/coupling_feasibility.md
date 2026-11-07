# DHSVM–MODFLOW coupling: feasibility and development assessment

Assessment date: 7 November 2026. Compute target: workstation. This is a
literature/source assessment and proposed experimental design, not a measured
runtime benchmark or a calibrated model. Hardware capacity and the pilot
watershed have not yet been specified.

## Feasibility by scale

An hourly coupled hydrologic analysis is a credible watershed-scale research
target. Rainier-scale operation is conditional on domain choice, groundwater
conceptualization, selective refinement and ensemble size. A coarse western
Washington analysis may be achievable on a workstation; a uniformly fine,
fully coupled regional ensemble cannot presently be promised. Core-code
optimization can improve throughput, but cannot supply missing observations or
make unresolved processes physically valid.

The following resolutions are starting experiments, not validated requirements
or guarantees. Refine vertically around soil horizons, confining layers and
moving water tables according to the application.

| Scale | Proposed initial configuration | Workstation assessment |
| --- | --- | --- |
| Pilot watershed, approximately 100–1,000 km² | DHSVM at 50–150 m; groundwater at 100–500 m with local refinement and approximately 5–15 hydrostratigraphic layers; hourly exchange with event substeps | Most credible first target. Begin with one member, then test 8, 16 and 32 members. |
| Rainier and selected surrounding drainages | Surface at 100–250 m; groundwater at 250–1,000 m, with finer valley/hillslope submodels and a modest ensemble | Plausible after pilot benchmarking. Groundwater structure and glacier routing may limit scientific accuracy more than surface-model runtime. |
| Western Washington | Surface at 500–1,000 m; multiple connected hydrogeologic domains at approximately 500–2,000 m, with detailed nests for priority areas | Coarse hourly nowcasts merit testing. Uniform fine-resolution regional ensembles and multi-day forecasts issued hourly would be substantially more demanding. |

Rainier does not define one watershed. Select outlets in the Carbon, White,
Cowlitz, Nisqually and/or Puyallup systems as appropriate. Define groundwater
boundaries from geology and hydraulic connectivity, including possible
cross-divide flow, rather than automatically applying the surface mask.
[NPS drainage context](https://www.nps.gov/mora/planyourvisit/park-brochure.htm).

For an illustrative scaling calculation only, a 500 km² pilot, 3,000 km²
multi-basin Rainier domain and 70,000 km² regional domain would contain:

| Assumed domain area | Surface cells at 100 m | Surface cells at 500 m | Groundwater cells at 500 m and 10 layers |
| --- | ---: | ---: | ---: |
| 500 km² | 50,000 | 2,000 | 20,000 |
| 3,000 km² | 300,000 | 12,000 | 120,000 |
| 70,000 km² | 7,000,000 | 280,000 | 2,800,000 |

These are `area / cell_area` calculations, not delineated domain sizes. They
exclude inactive cells, nesting and boundary extensions. With an assumed 32
members, 20 scalar surface states per cell and float64 storage, seven million
surface cells alone require about 35.8 GB for ensemble states. At 500 m the
same illustrative allocation is about 1.43 GB. Neither figure includes solver
arrays, sparse matrices, groundwater states, parameter ensembles, forcing or
assimilation copies; these are not measured DHSVM memory footprints.

Parallel DHSVM has been demonstrated at 90 m over 25,000 and 668,000 km², with
maximum reported speedup of 105 on 480 processors. This establishes that the
surface model can scale on distributed-memory machines; it does not benchmark
the proposed coupled ensemble on a workstation.
[Perkins et al., 2019](https://doi.org/10.1016/j.envsoft.2019.104533).

## Physical architecture

DHSVM is attractive for mountain snow, canopy, energy and hillslope hydrology.
MODFLOW supplies a more general multilayer aquifer representation. DHSVM
already simulates soil moisture and shallow saturated flow, however, so the
coupling must explicitly replace or partition overlapping processes.
[PNNL DHSVM source and model description](https://github.com/pnnl/DHSVM-PNNL).

Two architectures should be compared:

1. **Separate shallow and deep domains.** Retain DHSVM shallow hillslope
   hydrology, transfer leakage through a distinct intervening vadose domain to
   MODFLOW, and return deeper groundwater feedback where material. This can be
   economical where a physical separation exists. It requires an explicit
   partition and is unsuitable if the two codes describe the same aquifer.
2. **One continuous saturated system.** Retain DHSVM snow/canopy/land-surface
   processes and the supported part of its soil calculation; replace its
   overlapping saturated routing with MODFLOW. Couple the vadose/soil solver
   to groundwater head, exfiltration and infiltration constraints. This is
   the more appropriate research direction where the water table enters the
   root zone, but requires a more substantial DHSVM refactor.

In either design, assign stream routing/exchange and ET once. A possible
implementation routes land runoff/interflow into MODFLOW SFR, with stream
gain/loss solved against groundwater. Other routing choices are possible,
provided their exchanges and storage are conservative. A land-surface component
must not remove water for ET that a second component removes again.

Use a separate vadose module only for physical storage not already represented
by DHSVM. Its top boundary is the actual interface drainage, not automatically
land-surface rainfall. Return upward feedback where the selected equations
support it; an imposed downward recharge flux cannot reproduce all shallow
water-table/soil interactions.

The main physics choice is the vadose approximation. Standard MODFLOW 6 UZF
uses gravity-dominated vertical kinematic flow rather than full Richards flow.
The USGS manual warns about neglected capillary gradients for short intense
infiltration events and omits an explicit capillary fringe. Benchmark its
adequacy for the intended wetting fronts and shallow groundwater; hazard
studies involving suction or perched pressures may need Richards-based local
models. [USGS groundwater-flow theory, UZF limitations](https://pubs.usgs.gov/tm/06/a55/tm6a55.pdf).

An October 2026 Richards implementation for MODFLOW 6 is a candidate benchmark,
but its authors identify it as research code outside official USGS releases.
Evaluate its physics, maintenance and performance before adopting it.
[Paper](https://doi.org/10.1111/gwat.70112),
[author implementation and benchmarks](https://www.gsienv.com/software/modflow-6-richards/modflow6-richards/).

## Existing precedents and alternatives

- Lowry et al. coupled DHSVM to a local COMSOL groundwater model, demonstrating
  component feasibility for mountain meadow groundwater. It is not a maintained
  DHSVM–MODFLOW operational system.
  [Lowry et al., 2011](https://doi.org/10.1029/2010WR010086).
- A 2026 DHSVM extension adds dynamic bidirectional stream–shallow-groundwater
  exchange. It is a useful comparison baseline, while retaining a medium
  complexity groundwater approximation rather than a deep 3-D aquifer solver.
  [Boardman et al., 2026](https://doi.org/10.5194/hess-30-5145-2026).
- GSFLOW is a well-established PRMS–MODFLOW integration precedent, but its
  official implementation operates daily and therefore does not directly meet
  this hourly objective.
  [USGS GSFLOW description](https://www.usgs.gov/software/gsflow-coupled-groundwater-and-surface-water-flow-model).
- Current WRF-Hydro documentation includes a two-way MODFLOW 6 coupling example
  exchanging bottom drainage and groundwater drain flux. It is a useful
  engineering comparator if surface-model choice remains open; this example
  alone does not establish adequate Richards physics or hazard validation.
  [WRF-Hydro appendix A20.7](https://wrf-hydro.readthedocs.io/en/latest/appendices.html#a20-7-example-coupling-wrf-hydro-modflow-6).

Targeted searches did not identify a maintained turnkey DHSVM–MODFLOW 6
coupler. The proposed system should be budgeted as new coupling research and
software development, with existing interfaces reused wherever possible.

## Rainier and regional data constraints

The Upper White River is a useful candidate for assessing glacier–stream–
groundwater interactions: a USGS study provides seepage measurements and a
groundwater interpretation. Its limitations explicitly identify inaccurate
glacier-melt recharge in the underlying SWB approach. Snow-only recharge would
therefore be an inadequate default for a glacier-influenced pilot.
[USGS Upper White River study](https://pubs.usgs.gov/publication/sir20245015/full).

The existing rainier3d work provides relevant GIS/source infrastructure and
historical wells. Its paper explicitly states that it does not solve groundwater
flow or assimilate hydrologic observations. Static well levels and global
water-table estimates are useful priors, not an hourly monitoring network.
[rainier3d primary project paper](https://denolle-lab.github.io/mt-rainier-digital-model/paper/).

For the lowlands, the USGS southeastern Puget Sound model provides a calibrated
13-layer, 500-ft MODFLOW-NWT framework with monthly 2005–2015 transients. Reuse
its hydrostratigraphy and historical constraints where applicable, while
revalidating an MF6 conversion and finer temporal discretization. It is neither
an all-western-Washington model nor demonstrated hourly forecasting software.
[USGS numerical model](https://doi.org/10.3133/sir20245026v2).

Western Washington includes distinct volcanic, fractured-bedrock, glacial and
alluvial systems. A regional architecture needs appropriate geological
provinces and connected boundaries, with local models for detailed hazard
applications. [USGS Puget aquifer framework](https://doi.org/10.3133/pp1424D).

## Assimilation and near-real-time operation

Begin with offline spin-up/calibration, then add a localized ensemble filter
such as an ensemble square-root filter or LETKF. Benchmark ensemble size rather
than assuming 32 members are sufficient. Surface moisture observations require
depth/support-aware operators; wells require screen- and datum-aware head
operators; discharge needs upstream routing and storage. Treat precipitation
as uncertain forcing and distinguish shared errors when gauges also enter
gridded forcing products. [UCAR ensemble DA principles](https://docs.dart.ucar.edu/en/stable/guide/introduction-ensemble-da.html),
[PDAF supported filters](https://www.awi.de/en/science/special-groups/scientific-computing/data-assimilation/pdaf-parallel-data-assimilation-framework.html).

The proposed observations complement one another, but do not uniquely constrain
the deep subsurface. Add snow/SWE constraints to separate snow accumulation and
melt errors from recharge errors. Prioritize nested piezometers, soil moisture
profiles and pressure/suction measurements near hazard-relevant depths when
designing new instrumentation. Temperature, humidity, wind and radiation or
defensible estimates are also needed for snow energy calculations.

Track every assimilation-induced storage increment and reconcile moisture and
pressure states physically. Improved surface fit can otherwise accompany worse
deep-soil estimates or inconsistent fluxes.
[Li et al., mass-constrained assimilation](https://doi.org/10.5194/hess-16-105-2012).

Hourly output, hourly analysis and low observation latency are separate goals:

- USGS continuous stations may provide frequent telemetry, but inventory actual
  site availability and latency for streamflow and groundwater separately.
  [USGS latest-continuous API](https://api.waterdata.usgs.gov/ogcapi/v0/collections/latest-continuous?f=html).
- Ecology flow gages can log at 15/30 minutes while transmitting every three
  hours. [Ecology DataStream](https://apps.ecology.wa.gov/continuousflowandwq/).
- NOAA's current hourly multisensor MRMS table lists 20-minute latency for
  Pass 1 and one-hour latency for Pass 2. Actual ingestion adds further delay.
  [NOAA MRMS operational product table](https://www.nssl.noaa.gov/projects/mrms/operational/tables.php).
- Paradise SNOTEL supplies hourly snow/weather variables, but sensor sets and
  reporting latency must be checked station by station; do not assume an hourly
  soil-moisture network. [NRCS Paradise station](https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=679&state=).

Publish an hourly provisional nowcast using available observations and preserve
a rolling restart window for late-data replay. Report valid time and issue
time. For low-latency current-time estimates, forecast from the most recent
available forcing/analysis and label that interval as a forecast. Multi-day
forecasts reissued hourly require a separate compute budget.

## Core-code and memory optimization

Use Python initially for scheduling and configuration, with compiled model
kernels and sparse/vectorized conservative grid mapping. Optimize measured
bottlenecks, preserving the coupled equations and hazard-relevant accuracy.

| Order | Action | Implementation and qualification |
| --- | --- | --- |
| 1 | Instrument and reduce unnecessary I/O | MF6 already supports profiling/memory reporting and selective head/budget output. Measure wall-clock time externally as well. |
| 2 | Keep models alive and exchange arrays in memory | Use the MF6 API. Refactor DHSVM initialization, advance, state/flux access, checkpoint and finalize into an explicit context interface. |
| 3 | Remove repeated allocation and copies | Inspected DHSVM `RouteSubSurface.c` allocates routing work arrays each call. Test reusable contiguous buffers and configurable diagnostic output. Cache only invariant geometry. |
| 4 | Schedule ensembles to fit actual RAM/cores | Compare independent member processes with within-model parallelism; avoid oversubscribing numerical libraries. Localize the filter and batch work. |
| 5 | Tune and optimize numerical kernels | Profile solver/preconditioner choices, array layout and vectorization before changing algorithms. Parallel routing needs race-free neighbor/channel accumulation. |
| 6 | Evaluate advanced alternatives | Adaptive deep-aquifer updates, reduced precision, GPU kernels and surrogates require measured gains and error tests; none is a default guarantee. |

The MODFLOW API already permits fine-grained coupling; rewriting a second API
is unnecessary. Its documented state/pointer and lifecycle semantics must be
verified for the selected release. Use isolated persistent workers when shared
module state prevents safe concurrent instances. Also verify actual time-step
control rather than assuming an API argument overrides TDIS/ATS.
[API paper](https://doi.org/10.1016/j.envsoft.2021.105257),
[BMI source reference](https://modflow-org.github.io/modflow6/namespacemf6bmi.html),
[XMI source reference](https://modflow-org.github.io/modflow6/namespacemf6xmi.html).

The repeated DHSVM allocations and diagnostic writes are visible in the public
source, but their runtime contribution has not been measured for this project.
Do not attach a speedup estimate to these changes yet.
[Inspected routing source](https://raw.githubusercontent.com/pnnl/DHSVM-PNNL/master/DHSVM/sourcecode/RouteSubSurface.c),
[MF6 profiling controls](https://modflow6.readthedocs.io/en/latest/_mf6io/sim-nam.html),
[IMS solver controls](https://modflow6.readthedocs.io/en/latest/_mf6io/sln-ims.html).

## Hazard applications

An hourly hydrologic backbone can supply evolving initial conditions for all
requested hazards. The actual hazard calculations require distinct physics
and resolutions, as specified in the project's `AGENTS.md`.

- Earthquake liquefaction/lateral spreading need shaking, susceptible
  stratigraphy and strength/resistance or deformation models. Hourly
  groundwater conditions initialize the event; they do not resolve cyclic
  pressure generation or displacement.
  [USGS earthquake effects](https://www.usgs.gov/programs/earthquake-hazards/what-are-effects-earthquakes).
- Landslide analyses need pressure/suction at relevant depths plus stability
  mechanics. Compare representative hillslopes against an appropriate
  infiltration/stability reference before regional extrapolation.
  [USGS TRIGRS](https://pubs.usgs.gov/of/2008/1159/).
- Flood discharge is a forcing/input for inundation hydraulics. Event-scale
  peaks and flash floods can require subhourly forcing and routing; fine DEM
  sampling does not create fine-scale hydraulic validation.
- Wildfire applications require a separate fuel/fire interpretation. Postfire
  runoff and debris-flow initiation can depend on 5–30-minute rainfall bursts;
  retain that information rather than averaging it away.
  [USGS postfire analysis](https://pubs.usgs.gov/publication/ofr20231025/full).
- Burn impacts on infiltration are not universal. A recent western Oregon
  Cascades study found higher field-saturated conductivity in burned soils,
  illustrating the need for regional evidence and competing scenarios.
  [Pimont et al., 2026](https://doi.org/10.1029/2025JG009611).

## Proposed development sequence and decision gates

1. **Select a monitored pilot and audit the workstation.** Define the outputs,
   basin/outlets, groundwater domain, geology, data latency and relevant hazard
   sites. Establish RAM/core/storage limits and a source/version inventory.
2. **Build physical references.** Validate snow/soil behavior and groundwater
   independently, then test a conservative one-way drainage coupling. Include
   vadose delay. Establish budgets and restart equivalence.
3. **Implement feedback.** Choose one of the explicit ownership architectures;
   test moving water tables, perched layers and gaining/losing streams against
   observations and suitable Richards-based local cases.
4. **Profile and optimize.** Test multiple grids and coupling intervals, first
   one member then 8/16/32. Measure peak RAM, I/O, convergence, failed steps and
   p50/p95/p99 end-to-end cycle times during demanding events. Keep a pinned
   scientific reference for every optimization.
5. **Add assimilation and replay.** Compare open-loop and assimilated forecasts
   at withheld sites/events. Require physically consistent increments and
   demonstrable forecast improvement, not only better analysis fit.
6. **Expand conditionally.** Move to additional Rainier basins after pilot
   accuracy and latency targets pass. Build a coarse regional backbone with
   detailed nests; reassess hardware if the required detail exceeds workstation
   capacity.

A reasonable proposed engineering target is completing recurring nowcast
processing within 15–30 minutes, leaving recovery/replay margin in an hourly
cycle. This is an acceptance target, not a runtime forecast. Data arrival delay
adds to product freshness even if computational throughput meets the target.
Benchmark longer forecasts independently: under similar per-step cost, a
72-hour forecast requires approximately 72 times the propagation work of a
one-hour nowcast, before accounting for different I/O and assimilation work.

The current repository's monthly PRISM/recharge and MF6 template provide
configuration/data-management scaffolding. Implementing this architecture
requires new hourly meteorological ingestion, physical coupling, groundwater
setup, full-state persistence, assimilation and performance instrumentation.
No coupled runtime or hazard skill has been demonstrated by the present
assessment.
