"""A synthetic 3D staircase valley exercising the Richards reference kernel.

Uniform, deliberately weak rainfall is imposed as infiltration on horizontal
top faces. Every accepted step checks the pressure needed by that flux: the
experiment fails if a ponded boundary would be needed. There is no runoff,
channel, snow, ET or observation assimilation in this subsurface experiment.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import time

import numpy as np

from hydrostate.flow import BoundaryFaces, Grid3D, RichardsModel, StepResult


@dataclass(frozen=True)
class ValleyConfig:
    nx: int = 10
    ny: int = 8
    nz: int = 12
    length_x_m: float = 100.0
    length_y_m: float = 80.0
    height_m: float = 6.0
    hours: int = 72
    max_step_s: float = 600.0

    def __post_init__(self) -> None:
        if min(self.nx, self.ny, self.nz) < 2 or self.hours < 1:
            raise ValueError("At least two cells on each axis and a positive duration are required")
        if min(self.length_x_m, self.length_y_m, self.height_m, self.max_step_s) <= 0:
            raise ValueError("Domain dimensions and time step must be positive")


def build_valley(config: ValleyConfig | None = None) -> tuple[RichardsModel, np.ndarray, np.ndarray, np.ndarray]:
    """Model, initial pressure, actual stair surface [y,x], and material labels.

    The target valley floor declines eastward and rises on both sides. Cells
    are entirely retained below that surface, producing real Cartesian stair
    faces. No terrain-following derivatives or rotated-coordinate shortcuts.
    """
    c = config or ValleyConfig()
    dx, dy, dz = c.length_x_m / c.nx, c.length_y_m / c.ny, c.height_m / c.nz
    x, y, z = (np.arange(c.nx) + .5) * dx, (np.arange(c.ny) + .5) * dy, (np.arange(c.nz) + .5) * dz
    yy, xx = np.meshgrid(y, x, indexing="ij")
    target = 3.5 + .015 * (c.length_x_m - xx) + .025 * abs(yy - c.length_y_m / 2)
    active = z[:, None, None] + dz / 2 <= target[None, :, :]
    grid = Grid3D(x, y, z, dx, dy, dz, active)
    surface = np.max(np.where(active, grid.elevation_m + dz / 2, -np.inf), axis=0)
    depth = surface[None, :, :] - grid.elevation_m
    ks = np.full(grid.shape, 2e-5)
    ks[grid.elevation_m < 1.5] = 3e-6
    alpha = np.full(grid.shape, 1.2)
    theta_s = np.full(grid.shape, .43)
    n = np.full(grid.shape, 1.6)
    layer = np.full(grid.shape, "sandy_soil", dtype="<U24")
    layer[grid.elevation_m < 1.5] = "weathered_substrate"
    lens = (
        (depth >= 1.0) & (depth <= 1.5)
        & (xx[None, :, :] >= .3 * c.length_x_m) & (xx[None, :, :] <= .7 * c.length_x_m)
        & (yy[None, :, :] >= .25 * c.length_y_m) & (yy[None, :, :] <= .75 * c.length_y_m)
    )
    ks[lens] = 3e-7
    alpha[lens] = .6
    theta_s[lens] = .46
    n[lens] = 1.4
    layer[lens] = "low_K_lens"
    layer[~active] = "inactive"
    model = RichardsModel(grid, ks_m_s=ks, alpha_per_m=alpha, theta_s=theta_s, n=n)
    initial_head = 2.3 + .003 * (c.length_x_m - xx) + .006 * abs(yy - c.length_y_m / 2)
    psi = np.where(active, initial_head[None, :, :] - grid.elevation_m, np.nan)
    return model, psi, surface, layer


def rainfall_mm_h(hour: int) -> float:
    """Hypothetical uniform vertical liquid rain, not measured forcing."""
    return .05 + (.15 if 12 <= hour < 30 else 0.0) + (.1 if 48 <= hour < 60 else 0.0)


def valley_boundaries(model: RichardsModel, rain_mm_h: float) -> tuple[BoundaryFaces, BoundaryFaces]:
    """Regional head only on the outer east plane; internal stair risers sealed.

    Grid3D.boundary_cells deliberately returns ALL exposed faces. An internal
    east-facing topographic cliff is not the regional model outlet.
    """
    g = model.grid
    east = g.boundary_cells(2, 1)
    east = east[np.unravel_index(east, g.shape)[2] == g.shape[2] - 1]
    top = g.boundary_cells(0, 1)
    density = rain_mm_h / 1000 / 3600
    return (
        BoundaryFaces(east, 2, 1, "head", 2.3),
        BoundaryFaces(top, 0, 1, "flux", -density * g.face_area_m2(0)),
    )


def _json_field(array: np.ndarray) -> list:
    """Use JSON null outside the physical domain, never nonstandard NaN."""
    return np.where(np.isfinite(array), array, None).tolist()


def run_valley(
    *, config: ValleyConfig | None = None, hours: int | None = None,
    snapshot_interval_h: int = 1,
) -> dict:
    """Run a nonponding 3D subsurface demonstration and return JSON-ready data."""
    settings = asdict(config or ValleyConfig())
    if hours is not None:
        settings["hours"] = hours
    c = ValleyConfig(**settings)
    if snapshot_interval_h < 1:
        raise ValueError("snapshot_interval_h must be positive")
    model, psi, surface, layer = build_valley(c)
    g = model.grid
    top = g.boundary_cells(0, 1)
    initial_storage = float(np.nansum(model.storage(psi)) * g.cell_volume_m3)
    area = c.length_x_m * c.length_y_m
    precipitation = outlet = 0.0
    records, snapshots, cycle_times = [], [], []
    accepted = rejected = iterations = 0
    max_residual = max_increment = max_x_flux = max_y_flux = 0.0
    max_top_pressure = -np.inf
    min_step = c.max_step_s
    start = datetime(2026, 1, 1, tzinfo=timezone.utc)
    started = time.perf_counter()

    def store(hour: int, rain: float, q: float, count: int, wall: float) -> None:
        water = float(np.nansum(model.storage(psi)) * g.cell_volume_m3)
        theta = model.theta(psi)
        top_theta = theta.ravel()[top]
        residual = water - initial_storage - precipitation + outlet
        records.append({
            "time_h": hour, "utc": (start + timedelta(hours=hour)).isoformat().replace("+00:00", "Z"),
            "precip_mm_h": rain, "outlet_m3_s": q,
            "mean_theta": float(np.nanmean(theta)), "theta_surface": float(top_theta.mean()),
            "soil_storage_m3": water, "budget_residual_m3": residual,
            "budget_residual_mm": residual / area * 1000,
            "substeps": count, "wall_time_s": wall,
        })
        if hour % snapshot_interval_h == 0 or hour == c.hours:
            snapshots.append({
                "time_h": hour, "theta": _json_field(theta),
                "pressure_head_m": _json_field(psi),
                "total_head_m": _json_field(psi + g.elevation_m),
            })

    store(0, 0.0, 0.0, 0, 0.0)
    for hour in range(c.hours):
        cycle_start = time.perf_counter()
        rain = rainfall_mm_h(hour)
        flux_density = rain / 1000 / 3600
        head_boundary, infiltration = valley_boundaries(model, rain)

        def audit(result: StepResult, end_s: float) -> None:
            nonlocal max_top_pressure, max_x_flux, max_y_flux
            p = result.pressure_head_m
            kval = model.conductivity(p).ravel()[top]
            # Reconstruct the pressure required by this prescribed Neumann
            # flux using the same half-cell Darcy conductance as the FV grid.
            implied = p.ravel()[top] - g.dz_m / 2 + flux_density * g.dz_m / (2 * kval)
            maximum = float(implied.max())
            max_top_pressure = max(max_top_pressure, maximum)
            if maximum > 1e-6:
                raise ValueError("Prescribed infiltration requires ponding; this nonponding experiment is invalid")
            ids = result.fluxes.internal_i
            js = result.fluxes.internal_j
            values = abs(result.fluxes.internal_m3_s)
            xface, yface = (js - ids) == 1, (js - ids) == c.nx
            max_x_flux = max(max_x_flux, float(values[xface].max(initial=0)))
            max_y_flux = max(max_y_flux, float(values[yface].max(initial=0)))

        result = model.advance_interval(
            psi, 3600, boundaries=(head_boundary, infiltration), max_step_s=c.max_step_s,
            tolerance_m3_s=1e-10, head_tolerance_m=1e-7, on_accept=audit,
        )
        psi = result.pressure_head_m
        rain_volume = flux_density * area * 3600
        east_volume = result.boundary_out_volume_m3 + rain_volume
        precipitation += rain_volume
        outlet += east_volume
        accepted += result.accepted_steps
        rejected += result.rejected_steps
        iterations += result.nonlinear_iterations
        max_residual = max(max_residual, result.max_local_residual_m3_s)
        max_increment = max(max_increment, result.max_pressure_increment_m)
        min_step = min(min_step, result.min_step_s)
        cycle = time.perf_counter() - cycle_start
        cycle_times.append(cycle)
        store(hour + 1, rain, east_volume / 3600, result.accepted_steps, cycle)
    final_storage = float(np.nansum(model.storage(psi)) * g.cell_volume_m3)
    signed = final_storage - initial_storage - precipitation + outlet
    denominator = max(precipitation, abs(outlet), 1e-12)
    max_ledger = max(abs(row["budget_residual_m3"]) for row in records)
    return {
        "meta": {
            "schema_version": 1, "model": "Synthetic 3D mixed-Richards reference valley",
            "config": asdict(c), "start_utc": start.isoformat().replace("+00:00", "Z"),
            "coordinate_reference": {"horizontal_crs": "synthetic local Cartesian metres", "vertical_datum": "arbitrary local elevation datum"},
            "array_order": "[z,y,x], z cell centres ascending and positive upward",
            "aquifer_id": "synthetic_unconfined", "precipitation_spatial_support": "uniform",
            "units": "SI; pressure head m, total head=z+psi m, water content m3/m3, rainfall mm/hour, boundary discharge m3/s",
            "outlet": {"id": "east_valley_boundary", "kind": "fixed_total_head_boundary", "x_m": c.length_x_m, "y_m": c.length_y_m / 2, "total_head_m": 2.3},
            "record_support": "States at end of hour; rainfall and net east-boundary discharge averaged over preceding hour",
            "equation": "V*d[theta(psi)+Ss*max(psi,0)]/dt + sum(outward Darcy face-volume fluxes)=0; H=z+psi",
            "boundaries": [
                "Uniform liquid infiltration on horizontal top faces; required boundary suction checked at every accepted substep",
                "East vertical model boundary has fixed total hydraulic head 2.3 m on every active exposed face",
                "Bottom, west, north/south and vertical staircase faces are impermeable",
            ],
            "assumptions": [
                "Hypothetical Cartesian 3D valley and material values; all coordinates and dates are synthetic",
                "Active cells form a staircase valley surface; no cut cells, terrain metric approximation, or claimed mesh convergence",
                "Unified isotropic van Genuchten--Mualem Richards continuum; harmonic internal face K and small saturated elastic storage",
                "Sources are prescribed infiltrating liquid rain only: snow/glacier/frozen soil, canopy, roots/ET, surface ponding/runoff and channels are omitted",
                "Weak prescribed infiltration is accepted only while the discrete boundary-pressure reconstruction remains unponded",
                "East-boundary discharge is an open-groundwater-boundary surrogate, not a calibrated river-discharge prediction",
                "No observations assimilated; no earthquake, liquefaction, landslide, wildfire, flood or other hazard predictions",
                "No hysteresis, preferential/fracture flow, air dynamics or temperature-dependent transport",
            ],
            "parameter_provenance": "Illustrative values: soil Ks=2e-5m/s, substrate=3e-6m/s, lens=3e-7m/s; not measured or calibrated",
        },
        "grid": {
            "x_m": g.x_m.tolist(), "y_m": g.y_m.tolist(), "z_m": g.z_m.tolist(),
            "active": g.active.tolist(), "surface_elevation_m": surface.tolist(),
            "dx_m": g.dx_m, "dy_m": g.dy_m, "dz_m": g.dz_m,
            "layer_name": layer.tolist(),
        },
        "snapshots": snapshots, "records": records,
        "budget": {
            "initial_storage_m3": initial_storage, "final_storage_m3": final_storage,
            "precip_m3": precipitation, "outlet_net_m3": outlet,
            "signed_residual_m3": signed, "signed_residual_mm": signed / area * 1000,
            "max_abs_cumulative_residual_m3": max_ledger,
            "normalized_error": abs(signed) / denominator,
            "normalization_denominator_m3": denominator,
            "normalization_definition": "max(precipitation volume, absolute net east boundary volume); initial storage excluded",
            "equation": "final storage - initial storage - precipitation + net east boundary outflow = residual",
            "acceptance_tolerance_mm": .001,
            "accepted": max_ledger / area * 1000 < .001,
        },
        "diagnostics": {
            "wall_time_s": time.perf_counter() - started,
            "timing_scope": "forcing, solver, audits and in-memory output construction; excludes import and JSON writing",
            "accepted_substeps": accepted, "rejected_substeps": rejected,
            "nonlinear_iterations": iterations, "min_accepted_step_s": min_step,
            "max_local_residual_m3_s": max_residual, "local_tolerance_m3_s": 1e-10,
            "max_pressure_increment_m": max_increment, "pressure_tolerance_m": 1e-7,
            "max_implied_top_pressure_head_m": max_top_pressure,
            "max_abs_internal_x_flux_m3_s": max_x_flux,
            "max_abs_internal_y_flux_m3_s": max_y_flux,
            "hour_cycle_p95_s": float(np.percentile(cycle_times, 95)),
            "active_cells": int(g.active.sum()), "mesh_convergence_established": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hours", type=int, default=72)
    parser.add_argument("--output", type=Path, default=Path("runs/valley3d/results.json"))
    parser.add_argument("--snapshot-interval-h", type=int, default=1)
    args = parser.parse_args(argv)
    result = run_valley(hours=args.hours, snapshot_interval_h=args.snapshot_interval_h)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps({"output": str(args.output), "budget": result["budget"], "diagnostics": result["diagnostics"]}, indent=2))
    return 0 if result["budget"]["accepted"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
