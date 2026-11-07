# Prompt for Antigravity: executable integration and operational integrity

Copy the text below as a standalone review request.

---

Act as an independent adversarial numerical/software reviewer. Determine whether
these projects can actually execute a reproducible observation-to-soil-state-to-
hazard forecast, and build small counterexamples that break unsupported claims.
Review the code and equations, not just documentation or passing test counts.
Keep upstream repositories unchanged; provide proposed patches and regression tests
as review artifacts. Do not publish products or change cloud observational data.

Repositories and reviewed physical-code baselines:

- https://github.com/gaia-hazlab/gwl-space-time-smooth at
  `435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f`.
- https://github.com/gaia-hazlab/hydrostate at
  `c9328f331c6e6985054c0651bde15cc72235fb64`.
- https://github.com/gaia-hazlab/mechanicalstate at
  `b86d3b2e8a138f238e9db94dd7ece1c6c89db3eb`.
- Updated integration proposal:
  https://github.com/gaia-hazlab/hydrostate/blob/main/docs/soil_state_architecture.md
- The Hydrostate model-data extension audits gaia-cli commit
  `aa8ed3def5af4c32a4a1e601252fcee439344ac8`; verify actual current APIs and versions.

Record exact checkout commits and dirty states. Read repository instructions and
lockfiles. Separate the reviewed baselines from later changes. Use small synthetic
cases first; a dataset download or unavailable GPU must not prevent code-level
counterexamples. Core Hydrostate/Mechanicalstate computations support Python 3.11+;
the pinned Gaia/rioxarray interoperability stack requires Python 3.12+.

The intended system assimilates moisture probes, screened wells, river discharge,
seismic dv/v and GNSS-IR into a physically consistent, uncertain soil state. It
forecasts storms/seasons and damage/recovery, then supplies Landlab landslide,
liquefaction/lateral-spreading and flood consumers at each consumer's own support.
A workstation is the initial target. Regional 90 m products, local 3D process
patches and site-scale hazards may coexist. The proposal reuses the regional
soil twin's observation/orchestration layer instead of duplicating it; this is a
hypothesis to test, not an architectural constraint you must endorse.

Trace one real executable path through ingestion, native-support observation
operator, state update, material reconciliation, checkpoint, forecast and an actual
hazard component. Identify where that path stops. Classify every link as executed,
unit-tested only, format-compatible only, documented target or absent.

Required adversarial experiments:

1. **Combined conservation.** Independently reproduce the one-cell closed capillary
   transfer: P=PET=0, theta_wp=.10, theta_fc=.28, theta_sat=.45, root depth=1 m,
   Sy=.30, initial DTW=.50 m, init='wp', dt_days=1. Verify soil plus aquifer storage
   against external inputs/outputs, not just the soil ledger. Also test DTW=1.50 m, below the root-zone base. Inspect every internal
   transfer and first step. The existing audit diagnostic is in
   `hydrostate/scripts/review_gwl_capillary.py`; independently verify its interpretation.
   Test bounds, negative inputs, overflow and inactive cells without masking defects.

2. **True restart and cycling.** Compare uninterrupted and segmented runs containing
   initial snowpack, accumulation/melt, wetting/drying reversals, an off-hour earthquake,
   multiple recovery times, changing hydraulic materials, a sensor analysis and delayed
   data replay. Check all stores, pressures, modes, branches, reference frames, future
   events, member identities and cumulative water ledgers. Detect forecasts that reset
   soil moisture, groundwater anomaly, SWE or seismic/hysteresis reference state.
   Test both an existing numerical boundary and an off-grid requested stop; distinguish
   legitimate discretization changes from lost restart memory.

3. **Typed physical contracts.** Attempt DTW/head sign swaps, NAVD88/ellipsoid mixing,
   feet/metres, confined head masquerading as phreatic GWL, equal shapes with shifted
   coordinates or reversed y, missing CRS, incompatible masks and alpha in 1/kPa
   copied as 1/m. Test physical theta versus total hydraulic storage, theta>porosity
   hidden by clipping, and explicit material-only versus observation-induced water
   changes. A valid serializer must still reject an owner-state mismatch.

4. **Time and provenance.** Test mm per 6-hour/day/month interval versus mm/day and
   SI rates using exact bounds, leap months and irregular intervals. Test naive UTC,
   issue/valid/arrival confusion, stale observations, overlapping seismic windows,
   reference changes, duplicate events and future leakage from full-window well means
   or calibration. Do not average instantaneous jumps into interval means.

5. **Uncertainty integrity.** Follow the actual matrices/arrays. Determine whether
   BLUE and seismic inversion retain correlated observation error and whether the
   returned uncertainty is marginal variance, covariance access or coherent samples.
   Permute/duplicate/mismatch member IDs across water, damage, materials and forcing;
   the handoff must reject this. Test double-counting original measurements and
   already assimilated retrievals. Do not manufacture a hazard ensemble by independently
   sampling marginal sigma rasters or sorting members separately by variable.

6. **Spatial/consumer coupling.** Verify native versus delivered support and coarse
   driver labels. Test conservative recharge/runoff transfer across unequal and partial
   cells, watershed boundaries and masks; bilinear interpolation is not a water-volume
   proof. Instantiate the named Landlab or other consumer with a documented version;
   prove the supplied groundwater/recharge fields enter its equations. Parsing ASCII
   alone is insufficient. Test whether a changing DTW actually changes the exact
   liquefaction implementation. Trace runoff to routed discharge and hydraulic
   inundation rather than renaming a source term or boundary flux.

7. **Numerical adequacy.** Run meaningful hydraulic and material-exchange time refinement
   and spatial/vertical refinement on layered/perched/water-table cases. Check local
   as well as global budgets, field errors, flux timing and hazard-threshold errors.
   Distinguish exact exponential recovery from the order of the combined split solver.
   Test representable states, failed-step rollback and geometry-changing event rejection.

8. **Operational cost.** Profile a complete small ensemble cycle: ingestion, observation
   operators, analysis, forecast, checkpoint, delayed replay, mapping and hazard call.
   Measure wall-time distribution and peak resident memory, solver/coupling iterations
   and I/O. Inspect dense n_cell×n_obs intermediates, full-domain array copies, output
   volume, process safety and oversubscription. Give explicitly labeled memory/cost
   bounds for watershed, Rainier and western-Washington scenarios. Do not extrapolate
   a tiny single-member run into demonstrated regional real-time performance.

Deliver:

- A reproducibility manifest: commits, environment, commands, executed/skipped/failed
  tests and unavailable dependencies. Separate an existing passing suite from your
  newly exposed counterexamples.
- Ranked P0/P1/P2 findings with file/line, minimal runnable reproducer, expected versus
  actual values, scientific consequence and proposed repair. State unconfirmed risks
  separately and report concerns that testing disproved.
- A compatibility matrix for quantities, units, grid/reference, supports, clocks,
  stochastic objects, complete state ownership and real consumer APIs.
- Minimal isolated patches/regression tests for confirmed defects; do not weaken
  tolerances or change reference physics merely to make a test pass. Identify any
  correction that invalidates prior calibration and requires renewed validation.
- A smallest executable vertical slice and interface specification: reduced regional
  product versus 3D restart, forecast/analyze/checkpoint, member-aware hazard handoff,
  and readiness rejection when required physics or support are missing.
- A measured performance report with bounded projections and a release gate separating
  experimental soil-state products, screening, conditional triggering and operational
  hazard forecasts.

End with a direct verdict: what should be reused, what should be consolidated,
what cannot currently interoperate, and which three fixes unblock the most reliable
end-to-end demonstration. Prefer concrete evidence over architectural expansion.
