# Prompt for Fable: scientific architecture and hazard adequacy

Copy the text below as a standalone review request.

---

Act as an independent, adversarial scientific reviewer of a proposed GAIA soil-state
system. Your job is to find reasons it would give a physically wrong, statistically
unidentifiable or falsely confident hazard forecast, then recommend the smallest
scientifically defensible architecture. Do not defend the proposed design. Treat
README claims, roadmap status, previous reviews and this prompt as hypotheses,
not evidence. Distinguish demonstrated defects, conditional risks and untested work.

Read actual source, tests, technical reports and current primary literature in:

- https://github.com/gaia-hazlab/gwl-space-time-smooth — now titled GAIA Digital
  Twin of Soil; reviewed baseline `435445f17b90e1cdc4ee4e12a9e1cb17338f7e3f`.
- https://github.com/gaia-hazlab/hydrostate — physical-code review baseline
  `c9328f331c6e6985054c0651bde15cc72235fb64`.
- https://github.com/gaia-hazlab/mechanicalstate — physical-code review baseline
  `b86d3b2e8a138f238e9db94dd7ece1c6c89db3eb`.
- The updated proposal:
  https://github.com/gaia-hazlab/hydrostate/blob/main/docs/soil_state_architecture.md

Record the exact revisions you inspect, identify changes since these baselines,
and separate code capability from stated intent. The groundwater repository has
advanced beyond its older monthly interpolation description; inspect its daily
water budget, footprint BLUE/Matérn analysis, snow/hysteresis and forecast code.

The objective is a near-real-time and forecasted soil-state distribution over
watersheds, Mt Rainier and western Washington, initially using a workstation.
Inputs include precipitation/weather, depth-specific moisture probes, screened
wells, gauged discharge, seismic dv/v and GNSS-IR. The state should support Landlab
landslide applications, geotechnical liquefaction/lateral-spreading applications,
flood routing/inundation and wildfire/postfire effects. Hydrological evolution
spans storms to seasons; damage can jump during shaking or disturbance and recover
over days to years. Output cadence, internal integration step, native observation
support, forecast issue time and useful warning lead are separate requirements.

The candidate revision is to reuse the GAIA soil twin as regional observation/
analysis/orchestration host, use Hydrostate as a selectable conservative hydraulic
backend/reference and Mechanicalstate for material/petrophysical/recovery laws,
with Gaia as the data plane and named hazard consumers. Fine physical patches
would be introduced where needed; 90 m is a possible reporting ledger, not a
universal information scale. Challenge this allocation and compare it with
credible alternatives, including established integrated physical engines. Do not
assume a new uniform statewide 3D solver or another independent controller is needed.

Perform these investigations:

1. Reconstruct the governing balances and constitutive assumptions of each actual
   model. Check that root-zone, vadose, groundwater, snow, surface and channel
   stores each have one owner. Trace recharge, capillary rise, ET, runoff,
   interflow, baseflow and assimilation increments through the combined budget.
   Explicitly test a closed shallow-water-table column with capillary rise and
   zero rain/PET, including a table below the root-zone base. Check first-step
   accounting, parameter changes and clipping.
   Distinguish specific yield from retention drainage and elastic storage.

2. Test whether the observations can distinguish moisture, pressure, temperature,
   hydraulic hysteresis, rate-dependent capillarity and structural damage.
   Construct at least three observationally similar states with different hazard
   outcomes: perched versus regional groundwater; confined head versus shallow
   phreatic surface; and moisture-driven versus damage-driven dv/v. Explain which
   extra observations or priors resolve each case, and what remains unknowable.
   Audit seismic kernels/reference windows, GNSS-IR footprint/depth and nuisance
   effects, well screens, discharge integration and shared retrieval/calibration
   errors. Do not assimilate the same information twice through raw and derived data.

3. Examine how a posterior becomes a forecast. Require preserved soil/aquifer,
   snow, hysteresis, material, damage/recovery and forcing states. Check whether
   current forecasts actually continue an analysis or reinitialize. Freeze data,
   baselines, parameters and seismic references at issue time. Distinguish
   retrospective reanalysis from a valid operational hindcast.

4. Audit uncertainty. Does the implementation provide marginal variances, a full
   covariance representation or coherent member histories? Can the claimed hazard
   probabilities be identified from that object? Challenge independent pixel/time
   draws, shared-error independence, forcing-method spread called posterior spread,
   and dependence among water, damage, material strength and shaking. Specify the
   minimum valid probabilistic object for each consumer.

5. Build a hazard-specific requirements matrix. For each application state the
   predicted quantity, controlling soil variables and depth support, extra trigger/
   strength/geometry inputs, appropriate model, temporal/spatial/vertical scales,
   observable constraints and validation gate. Separate regional susceptibility,
   event-conditional triggering, displacement, runout and inundation. Treat proposed
   resolutions as hypotheses requiring convergence and independent evidence.

6. Verify the actual downstream equations/APIs. Does the selected Landlab component
   consume transient pore pressure or reconstruct wetness from recharge? Does the
   exact published Sanger/Geyin/Maurer map interface accept changing groundwater,
   or does it require a new validated coefficient/model pathway? How does runoff
   generation become routed reach discharge and inundation? Distinguish a readable
   raster from a physically executed coupling. Include thin susceptible layers,
   perched pressure, same-hourly-total/different-rainfall-burst cases, bank geometry
   and postfire changes of either conductivity sign.

7. Critique Hydrostate and Mechanicalstate as strongly as the regional model.
   Assess the limits of their synthetic grids, nonlinear space/time convergence,
   isotropic conductivity, prescribed stress/contact laws, diagnostic versus
   hydraulic memory, fixed porosity/geometry and event labels. Small-strain
   stiffness loss does not establish strength, liquefaction resistance or failure.
   A small water residual or fast toy run does not establish regional accuracy.

Deliver:

- A one-paragraph verdict: can this become the binding soil-state framework, and
  under exactly which restrictions? What should be reused, consolidated or replaced?
- The ten most consequential findings, ranked as P0 (invalidates a balance/claimed
  quantity), P1 (blocks a proposed use) or P2 (important limitation). For each give
  file/line or primary-source evidence, confidence, a falsification experiment and
  the smallest correction. State when a concern was disproved.
- The hazard requirements matrix and an explicit observability/nonuniqueness table.
- A revised ownership diagram and mathematical description of the minimal common
  state, assimilation transaction and uncertainty handoff.
- Five decisive acceptance experiments with quantitative observables and justified
  pass/fail criteria, including withheld-event forecast and probability calibration.
- A staged research plan identifying what must be fixed before any operational
  hazard claim and what can safely remain a regional screening approximation.

Use primary papers and official source/documentation. Report what you executed,
what you only inspected and what evidence was unavailable. Do not invent benchmark
results, treat calibration as independent validation or weaken a test to pass it.
