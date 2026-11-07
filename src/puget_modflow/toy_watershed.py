"""Synthetic, mass-conserving hillslope reference; not a DHSVM/MODFLOW coupler.

The finite-volume domain is a rotated rectangle. ``s`` points down a plane and
``d`` points into that plane; they are orthogonal physical distances, not an
untransformed terrain-following approximation. Elevation is
``z = z0 - s sin(beta) - d cos(beta)``. We solve the mixed Richards equation
``d[theta(psi) + Ss max(psi, 0)]/dt = div(K grad(z+psi)) - uptake``
with backward Euler and a conservative modified-Picard iteration. Hydraulic
conductivity follows van Genuchten--Mualem. The soil, surface detention and
channel reservoir are solved together. Fluxes are volume/time, positive out of
the indicated store. A small saturated compressibility permits positive heads.

All parameters and forcing are hypothetical. This is a single continuum,
two-dimensional reference experiment, not calibrated watershed prediction,
operational warning software, a glacier model, or a geomechanical model.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class ToyConfig:
    """SI geometry and numerical settings for an idealized headwater unit."""

    nx: int = 24
    nz: int = 24
    length_m: float = 120.0
    depth_m: float = 6.0
    width_m: float = 100.0
    slope_rise_run: float = 0.03
    surface_at_divide_m: float = 12.0
    hours: int = 120
    max_step_s: float = 600.0
    min_step_s: float = 0.5
    residual_tolerance_m3_s: float = 2.0e-9
    pressure_tolerance_m: float = 1.0e-5
    max_iterations: int = 70
    specific_storage_per_m: float = 1.0e-5
    initial_snow_mm: float = 25.0
    surface_detention_m: float = 0.001
    surface_tau_s: float = 1800.0
    channel_area_m2: float = 240.0
    channel_tau_s: float = 10800.0

    def __post_init__(self) -> None:
        if self.nx < 2 or self.nz < 4 or self.hours < 1:
            raise ValueError("nx >= 2, nz >= 4 and hours >= 1 are required")
        if min(self.length_m, self.depth_m, self.width_m, self.max_step_s) <= 0:
            raise ValueError("Geometry and time steps must be positive")
        if self.min_step_s <= 0 or self.min_step_s > self.max_step_s:
            raise ValueError("Require 0 < min_step_s <= max_step_s")
        if self.slope_rise_run < 0:
            raise ValueError("The toy's down-slope direction requires nonnegative slope")


def retention(
    pressure_head_m: np.ndarray,
    theta_r: np.ndarray | float = 0.06,
    theta_s: np.ndarray | float = 0.43,
    alpha_per_m: np.ndarray | float = 1.2,
    n: np.ndarray | float = 1.6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return theta, dtheta/dpsi and relative conductivity (Mualem l=0.5).

    ``psi < 0`` is matric suction; ``psi >= 0`` is saturated pressure. There is
    no hysteresis, ice, air-pressure dynamics, or preferential-flow domain.
    """
    psi = np.asarray(pressure_head_m, dtype=float)
    m = 1.0 - 1.0 / np.asarray(n)
    a = np.asarray(alpha_per_m) * np.maximum(-psi, 0.0)
    a_n = a ** n
    se = (1.0 + a_n) ** (-m)
    theta = np.asarray(theta_r) + (np.asarray(theta_s) - theta_r) * se
    capacity = (
        (np.asarray(theta_s) - theta_r)
        * m * n * np.asarray(alpha_per_m)
        * a ** (np.asarray(n) - 1.0)
        * (1.0 + a_n) ** (-m - 1.0)
    )
    # This algebra avoids cancellation in 1 - (1 - Se^(1/m))^m at dry states.
    bracket = 1.0 - (a_n / (1.0 + a_n)) ** m
    kr = np.sqrt(se) * bracket ** 2
    return theta, capacity, kr


class RichardsHillslope:
    """A conservative FV cross-section with optional top/channel boundaries.

    Left divide and bottom are impermeable. The top is a dry flux boundary or
    a ponded pressure boundary chosen by an active set. Each surface column
    has a linear detention reservoir above a threshold. Its runoff enters one
    channel reservoir (not a Saint-Venant channel/floodplain). At the toe,
    submerged faces exchange bidirectionally with that reservoir; exposed
    bank faces permit atmospheric seepage, never atmospheric water inflow.
    """

    def __init__(self, config: ToyConfig | None = None, *, lens: bool = True):
        self.config = config or ToyConfig()
        c = self.config
        self.dx = c.length_m / c.nx
        self.dz = c.depth_m / c.nz
        self.beta = float(np.arctan(c.slope_rise_run))
        self.cos_beta = float(np.cos(self.beta))
        self.sin_beta = float(np.sin(self.beta))
        self.x = (np.arange(c.nx) + 0.5) * self.dx
        self.depth = (np.arange(c.nz) + 0.5) * self.dz
        self.surface_z = c.surface_at_divide_m - self.x * self.sin_beta
        self.z = self.surface_z[None, :] - self.depth[:, None] * self.cos_beta
        self.volume = self.dx * self.dz * c.width_m
        self.surface_area = self.dx * c.width_m
        self.horizontal_area = c.length_m * c.width_m * self.cos_beta
        self.size = c.nx * c.nz
        self.theta_r = np.full((c.nz, c.nx), 0.06)
        self.theta_s = np.full((c.nz, c.nx), 0.43)
        self.alpha = np.full((c.nz, c.nx), 1.2)
        self.n = np.full((c.nz, c.nx), 1.6)
        self.ks = np.full((c.nz, c.nx), 1.5e-5)
        self.layer = np.full((c.nz, c.nx), "sandy_soil", dtype="<U24")
        weathered = np.broadcast_to(self.depth[:, None] >= 2.5, self.z.shape)
        self.ks[weathered] = 5.0e-6
        self.theta_s[weathered] = 0.38
        self.alpha[weathered] = 0.8
        self.layer[weathered] = "weathered_substrate"
        self.lens_mask = (
            (self.depth[:, None] >= 1.5) & (self.depth[:, None] < 2.25)
            & (self.x[None, :] >= 40.0) & (self.x[None, :] <= 90.0)
        ) if lens else np.zeros_like(self.z, dtype=bool)
        self.ks[self.lens_mask] = 6.0e-8
        self.alpha[self.lens_mask] = 0.45
        self.n[self.lens_mask] = 1.4
        self.theta_s[self.lens_mask] = 0.46
        self.layer[self.lens_mask] = "low_K_lens"
        # Root fractions are per model column and sum to exactly one.
        weights = np.exp(-self.depth / 0.40) * (self.depth <= 1.5)
        self.root_weights = weights / weights.sum()
        ids = np.arange(self.size).reshape(c.nz, c.nx)
        self.edge_i = np.r_[ids[:, :-1].ravel(), ids[:-1, :].ravel()]
        self.edge_j = np.r_[ids[:, 1:].ravel(), ids[1:, :].ravel()]
        self.edge_geom = np.r_[
            np.full(c.nz * (c.nx - 1), self.dz * c.width_m / self.dx),
            np.full((c.nz - 1) * c.nx, self.dx * c.width_m / self.dz),
        ]
        self.top_ids = ids[0, :]
        self.toe_ids = ids[:, -1]
        self.toe_z = (
            c.surface_at_divide_m - c.length_m * self.sin_beta
            - self.depth * self.cos_beta
        )
        self.initial_head = float(
            c.surface_at_divide_m - c.length_m * self.sin_beta - 1.8 * self.cos_beta
        )
        self.channel_bed = self.initial_head - 0.15
        self.sensor_iz = int(np.argmin(abs(self.depth - 1.25)))
        self.sensor_ix = int(np.argmin(abs(self.x - 65.0)))

    def theta(self, psi: np.ndarray) -> np.ndarray:
        return retention(psi, self.theta_r, self.theta_s, self.alpha, self.n)[0]

    def storage_density(self, psi: np.ndarray) -> np.ndarray:
        """Water volume per bulk volume, including small elastic storage."""
        return self.theta(psi) + self.config.specific_storage_per_m * np.maximum(psi, 0)

    def water_stress(self, psi: np.ndarray) -> np.ndarray:
        """Hypothetical Feddes-like dry and oxygen limits; not species calibrated."""
        dry = np.clip((psi + 60.0) / 57.0, 0.0, 1.0)
        oxygen = np.clip(-psi / 0.1, 0.0, 1.0)
        return dry * oxygen

    def initial_pressure(self, *, dry: bool = False) -> np.ndarray:
        pressure = self.initial_head - self.z
        if dry:
            # A nonequilibrium dry root zone over the same regional groundwater
            # body, not an unsupported dry channel or an altered channel bed.
            upper_unsaturated = (self.depth[:, None] <= 1.5) & (pressure < 0.0)
            pressure = np.where(upper_unsaturated, np.minimum(pressure, -10.0), pressure)
        return pressure

    def water_table_depth(self, psi: np.ndarray) -> list[float | None]:
        """Deepest base-connected zero-pressure crossing, excluding perches.

        Distances are normal to slope. If the base is unsaturated, the regional
        water table is unresolved below this finite domain and is returned null.
        """
        answer: list[float | None] = []
        for ix in range(self.config.nx):
            p = psi[:, ix]
            if p[-1] < 0:
                answer.append(None)
                continue
            iz = len(p) - 1
            while iz > 0 and p[iz - 1] >= 0:
                iz -= 1
            if iz == 0:
                answer.append(0.0)
            else:
                fraction = -p[iz - 1] / (p[iz] - p[iz - 1])
                answer.append(float(self.depth[iz - 1] + fraction * self.dz))
        return answer

    def _system(
        self, trial: np.ndarray, old: np.ndarray, dt: float,
        input_rate: float, pet_rate: float, *, closed: bool = False,
    ) -> tuple[Any, np.ndarray, dict[str, Any]]:
        """Residual and modified-Picard Jacobian in integrated m3/s units."""
        c = self.config
        count = self.size
        psi = trial[:count].reshape(self.z.shape)
        psi_old = old[:count].reshape(self.z.shape)
        surf = trial[count:count + c.nx]
        old_surf = old[count:count + c.nx]
        stage = trial[-1]
        theta, capacity, kr = retention(psi, self.theta_r, self.theta_s, self.alpha, self.n)
        density = theta + c.specific_storage_per_m * np.maximum(psi, 0.0)
        cap = capacity + c.specific_storage_per_m * (psi >= 0.0)
        head = (psi + self.z).ravel()
        kval = (self.ks * kr).ravel()
        harmonic = 2.0 * kval[self.edge_i] * kval[self.edge_j] / np.maximum(
            kval[self.edge_i] + kval[self.edge_j], 1.0e-30
        )
        conductance = self.edge_geom * harmonic
        qedge = conductance * (head[self.edge_i] - head[self.edge_j])
        residual = np.zeros(count + c.nx + 1)
        residual[:count] = (density - self.storage_density(psi_old)).ravel() * self.volume / dt
        np.add.at(residual, self.edge_i, qedge)
        np.add.at(residual, self.edge_j, -qedge)
        diagonal = np.zeros_like(residual)
        diagonal[:count] = np.maximum(cap.ravel(), 1.0e-12) * self.volume / dt
        np.add.at(diagonal, self.edge_i, conductance)
        np.add.at(diagonal, self.edge_j, conductance)
        rows = [self.edge_i, self.edge_j]
        cols = [self.edge_j, self.edge_i]
        data = [-conductance, -conductance]
        stress = self.water_stress(psi)
        uptake = (
            pet_rate * self.surface_area * self.root_weights[:, None] * stress
            if not closed else np.zeros_like(psi)
        )
        residual[:count] += uptake.ravel()
        qtop = np.zeros(c.nx)
        runoff = np.zeros(c.nx)
        qtoe = np.zeros(c.nz)
        if closed:
            residual[count:count + c.nx] = (surf - old_surf) * self.surface_area / dt
            diagonal[count:count + c.nx] = self.surface_area / dt
            residual[-1] = (stage - old[-1]) * c.channel_area_m2 / dt
            diagonal[-1] = c.channel_area_m2 / dt
        else:
            # Surface dry/ponded complementarity. At D=0 the exact finite-step
            # water supply is rain+melt plus any previous detention storage.
            gt = (
                0.5 * (kval[self.top_ids] + self.ks[0, :])
                * self.surface_area / (self.dz / 2.0)
            )
            available = self.surface_area * (input_rate + old_surf / dt)
            capacity_at_zero = gt * (self.surface_z - head[self.top_ids])
            ponded = (surf > 1.0e-10) | (available > capacity_at_zero + 1.0e-12)
            dry = ~ponded
            qtop[dry] = available[dry]
            qtop[ponded] = gt[ponded] * (
                self.surface_z[ponded] + surf[ponded] * self.cos_beta
                - head[self.top_ids[ponded]]
            )
            residual[self.top_ids] -= qtop
            residual[count:count + c.nx][dry] = surf[dry] * self.surface_area / dt
            diagonal[count:count + c.nx][dry] = self.surface_area / dt
            sid = count + np.arange(c.nx)[ponded]
            tid = self.top_ids[ponded]
            runoff[ponded] = (
                self.surface_area * np.maximum(surf[ponded] - c.surface_detention_m, 0.0)
                / c.surface_tau_s
            )
            runoff_deriv = (
                self.surface_area / c.surface_tau_s
                * (surf[ponded] > c.surface_detention_m)
            )
            residual[sid] = (
                self.surface_area * (surf[ponded] - old_surf[ponded]) / dt
                + qtop[ponded] + runoff[ponded] - input_rate * self.surface_area
            )
            diagonal[tid] += gt[ponded]
            diagonal[sid] += (
                self.surface_area / dt + gt[ponded] * self.cos_beta + runoff_deriv
            )
            rows.extend([tid, sid])
            cols.extend([sid, tid])
            data.extend([-gt[ponded] * self.cos_beta, -gt[ponded]])
            # The outlet is the actual bank ABOVE the reservoir bed. Faces
            # below the bed are no-flow, not a fully penetrating head boundary.
            # Split cut faces into submerged and exposed areas every iterate.
            gb = kval[self.toe_ids] * self.dz * c.width_m / (self.dx / 2.0)
            face_top = self.toe_z + self.dz / 2.0 * self.cos_beta
            face_bottom = self.toe_z - self.dz / 2.0 * self.cos_beta
            wetted_low = np.maximum(face_bottom, self.channel_bed)
            wetted_high = np.minimum(face_top, stage)
            wet_fraction = np.maximum(wetted_high - wetted_low, 0.0) / (self.dz * self.cos_beta)
            exposed_low = np.maximum(face_bottom, max(stage, self.channel_bed))
            exposed_height = np.maximum(face_top - exposed_low, 0.0)
            exposed_elevation = 0.5 * (face_top + exposed_low)
            exposed_fraction = exposed_height / (self.dz * self.cos_beta)
            wet_g = gb * wet_fraction
            exposed_g = gb * exposed_fraction * (head[self.toe_ids] > exposed_elevation)
            soil_g = wet_g + exposed_g
            qtoe = (
                wet_g * (head[self.toe_ids] - stage)
                + exposed_g * (head[self.toe_ids] - exposed_elevation)
            )
            residual[self.toe_ids] += qtoe
            diagonal[self.toe_ids] += soil_g
            rows.append(self.toe_ids)
            cols.append(np.full(c.nz, len(trial) - 1))
            data.append(-wet_g)
            volume = c.channel_area_m2 * (stage - self.channel_bed)
            old_volume = c.channel_area_m2 * (old[-1] - self.channel_bed)
            outlet = volume / c.channel_tau_s
            residual[-1] = (volume - old_volume) / dt + outlet - qtoe.sum() - runoff.sum()
            diagonal[-1] = (
                c.channel_area_m2 / dt + c.channel_area_m2 / c.channel_tau_s
                + wet_g.sum()
            )
            rows.extend([
                np.full(c.nz, len(trial) - 1),
                np.full(sid.size, len(trial) - 1),
            ])
            cols.extend([self.toe_ids, sid])
            data.extend([-soil_g, -runoff_deriv])
        rows.append(np.arange(len(trial)))
        cols.append(np.arange(len(trial)))
        data.append(diagonal)
        matrix = coo_matrix(
            (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
            shape=(len(trial), len(trial)),
        ).tocsr()
        diagnostics = {
            "uptake_m3_s": float(uptake.sum()), "qtop_m3_s": qtop,
            "qtoe_m3_s": float(qtoe.sum()), "runoff_m3_s": float(runoff.sum()),
            "local_max_residual_m3_s": float(np.max(np.abs(residual[:count]))),
        }
        return matrix, residual, diagnostics

    def advance(
        self, psi: np.ndarray, surface: np.ndarray, stage: float, dt: float,
        input_rate: float = 0.0, pet_rate: float = 0.0, *, closed: bool = False,
    ) -> tuple[np.ndarray, np.ndarray, float, dict[str, Any]]:
        """Advance one step or fail without modifying the beginning state."""
        c = self.config
        old = np.r_[psi.ravel(), surface, stage]
        trial = old.copy()
        last_pressure_increment = 0.0
        for iteration in range(1, c.max_iterations + 1):
            matrix, residual, diagnostics = self._system(
                trial, old, dt, input_rate, pet_rate, closed=closed
            )
            max_residual = float(np.max(np.abs(residual)))
            if (max_residual < c.residual_tolerance_m3_s
                    and last_pressure_increment < c.pressure_tolerance_m):
                if np.min(trial[self.size:self.size + c.nx]) < -1.0e-9:
                    raise RuntimeError("Nonphysical negative surface storage")
                if trial[-1] < self.channel_bed - 1.0e-9:
                    raise RuntimeError("Channel dried below its bed; reduce step/model needs dry-bed handling")
                diagnostics.update(
                    iterations=iteration, max_residual_m3_s=max_residual,
                    final_pressure_increment_m=last_pressure_increment,
                )
                return (
                    trial[:self.size].reshape(self.z.shape).copy(),
                    trial[self.size:self.size + c.nx].copy(),
                    float(trial[-1]), diagnostics,
                )
            increment = spsolve(matrix, -residual)
            if not np.all(np.isfinite(increment)):
                raise RuntimeError("Nonfinite Richards iterate")
            # Bound the *iteration* movement, not the accepted moisture or mass.
            # Iteration damping introduces no extra storage or flux.
            head_change = float(np.max(abs(increment[:self.size])))
            damping = min(1.0, 1.5 / max(head_change, 1.0e-30))
            trial += damping * increment
            last_pressure_increment = head_change * damping
        raise RuntimeError(f"Richards nonlinear solve did not converge in {c.max_iterations} iterations")


def synthetic_forcing(hour: int, scenario: str = "baseline") -> dict[str, float]:
    """Piecewise hourly *horizontal-area* forcing, independent of real dates."""
    if scenario not in {"baseline", "intense_storm", "postfire"}:
        raise ValueError(f"Unknown scenario: {scenario}")
    temperature = -2.0 if hour < 24 else 4.0 + 2.0 * np.sin(2 * np.pi * (hour - 8) / 24)
    precip = 0.8 if 6 <= hour < 18 else 0.0
    if 30 <= hour < 42:
        precip = 7.0
    if 60 <= hour < 66:
        precip = 12.0
    if 85 <= hour < 90:
        precip = 3.0
    if scenario == "intense_storm" and 60 <= hour < 66:
        precip *= 2.0
    daylight = max(0.0, float(np.sin(np.pi * ((hour % 24) - 6.0) / 12.0)))
    return {
        "precip_mm_h": precip,
        "temperature_c": float(temperature),
        "pet_mm_h": 0.18 * daylight,
    }


def run_toy(
    scenario: str = "baseline", *, config: ToyConfig | None = None,
    hours: int | None = None, max_step_s: float | None = None,
    snapshot_interval_h: int = 1,
) -> dict[str, Any]:
    """Run and return an entirely JSON-serializable synthetic experiment.

    Snapshot arrays use [normal depth][down-slope coordinate]. Every hourly
    record is an end-of-hour state; rates are preceding-hour means. A t=0
    record/snapshot is included. Snow and canopy are uniform over the slope.
    No observations, assimilation increments, DHSVM, or MODFLOW are invoked.
    """
    if scenario not in {"baseline", "intense_storm", "postfire"}:
        raise ValueError(f"Unknown scenario: {scenario}")
    if snapshot_interval_h < 1:
        raise ValueError("snapshot_interval_h must be positive")
    settings = asdict(config or ToyConfig())
    if hours is not None:
        settings["hours"] = hours
    if max_step_s is not None:
        settings["max_step_s"] = max_step_s
    c = ToyConfig(**settings)
    model = RichardsHillslope(c)
    if scenario == "postfire":
        # An explicitly hypothetical *one* postfire branch, not a universal law.
        # Reduced infiltration in the top 0.5 m and reduced interception/uptake.
        model.ks[model.depth <= 0.5, :] *= 0.3
    psi = model.initial_pressure()
    surface = np.zeros(c.nx)
    stage = model.initial_head
    snow = c.initial_snow_mm / 1000 * model.cos_beta
    canopy = 0.0
    canopy_capacity = (0.2 if scenario == "postfire" else 1.0) / 1000 * model.cos_beta
    slope_area = c.length_m * c.width_m
    def stores() -> dict[str, float]:
        return {
            "soil_m3": float(model.storage_density(psi).sum() * model.volume),
            "surface_m3": float(surface.sum() * model.surface_area),
            "channel_m3": float(c.channel_area_m2 * (stage - model.channel_bed)),
            "snow_m3": float(snow * slope_area),
            "canopy_m3": float(canopy * slope_area),
        }
    initial_stores = stores()
    initial_total = sum(initial_stores.values())
    accumulated = {"precip_m3": 0.0, "et_m3": 0.0, "outlet_m3": 0.0}
    records: list[dict[str, Any]] = []
    snapshots: list[dict[str, Any]] = []
    elapsed = 0.0
    rejected = 0
    total_substeps = 0
    max_iterations = 0
    max_local_residual = 0.0
    max_system_residual = 0.0
    max_accepted_pressure_increment = 0.0
    min_accepted_step = c.max_step_s
    max_hour_residual = 0.0
    start_utc = datetime(2026, 1, 1, tzinfo=timezone.utc)
    started = time.perf_counter()

    def save(hour: int, hourly: dict[str, float]) -> None:
        nonlocal max_hour_residual
        theta = model.theta(psi)
        table = model.water_table_depth(psi)
        root_theta = float(np.mean(np.sum(theta * model.root_weights[:, None], axis=0)))
        root_stress = float(np.mean(np.sum(model.water_stress(psi) * model.root_weights[:, None], axis=0)))
        storage = stores()
        residual = sum(storage.values()) - initial_total - accumulated["precip_m3"] + accumulated["et_m3"] + accumulated["outlet_m3"]
        max_hour_residual = max(max_hour_residual, abs(residual))
        record = {
            "time_h": hour, "utc": (start_utc + timedelta(hours=hour)).isoformat().replace("+00:00", "Z"),
            "precip_mm_h": hourly.get("precip_m3", 0.0) / model.horizontal_area * 1000,
            "rain_mm_h": hourly.get("rain_m3", 0.0) / model.horizontal_area * 1000,
            "melt_mm_h": hourly.get("melt_m3", 0.0) / model.horizontal_area * 1000,
            "snow_mm": snow / model.cos_beta * 1000,
            "temperature_c": hourly.get("temperature_c", synthetic_forcing(0, scenario)["temperature_c"]),
            "theta_root": root_theta, "wt_depth_m": table[-1],
            "groundwater_head_m": float((psi + model.z)[-1, -1]),
            "water_table_elevation_m": (
                float(model.surface_z[-1] - table[-1] * model.cos_beta) if table[-1] is not None else None
            ),
            "outlet_m3_s": hourly.get("outlet_m3", 0.0) / 3600,
            "root_stress": root_stress,
            "pressure_sensor_kpa": float(psi[model.sensor_iz, model.sensor_ix] * 9.80665),
            "budget_residual_m3": residual,
            "budget_residual_mm": residual / model.horizontal_area * 1000,
            "soil_storage_m3": storage["soil_m3"],
            "surface_storage_m3": storage["surface_m3"],
            "channel_storage_m3": storage["channel_m3"],
            "canopy_mm": canopy / model.cos_beta * 1000,
            "et_mm_h": hourly.get("et_m3", 0.0) / model.horizontal_area * 1000,
            "runoff_m3_s": hourly.get("runoff_m3", 0.0) / 3600,
            "groundwater_to_channel_m3_s": hourly.get("gw_channel_m3", 0.0) / 3600,
            "infiltration_mm_h": hourly.get("infiltration_m3", 0.0) / model.horizontal_area * 1000,
            "substeps": int(hourly.get("substeps", 0)),
            "solver_iterations": int(hourly.get("iterations", 0)),
        }
        records.append(record)
        if hour % snapshot_interval_h == 0 or hour == c.hours:
            snapshots.append({
                "time_h": hour, "theta": theta.tolist(), "pressure_head_m": psi.tolist(),
                "total_head_m": (psi + model.z).tolist(),
                "pore_pressure_pa": (psi * 9806.65).tolist(),
                "surface_water_depth_m": surface.tolist(),
                "water_table_depth_m": table,
                "channel_stage_m": stage,
            })

    save(0, {})
    next_step = c.max_step_s
    for hour in range(c.hours):
        forcing = synthetic_forcing(hour, scenario)
        hourly = {key: 0.0 for key in (
            "precip_m3", "rain_m3", "melt_m3", "et_m3", "outlet_m3",
            "runoff_m3", "gw_channel_m3", "infiltration_m3", "substeps", "iterations",
        )}
        hourly["temperature_c"] = forcing["temperature_c"]
        hour_end = (hour + 1) * 3600.0
        while elapsed < hour_end - 1.0e-8:
            dt = min(next_step, hour_end - elapsed)
            precip_rate = forcing["precip_mm_h"] / 1000 / 3600 * model.cos_beta
            rain_rate = precip_rate if forcing["temperature_c"] > 0 else 0.0
            snowfall_rate = precip_rate - rain_rate
            # Empirical 5 mm water / (deg C day), capped by available SWE.
            potential_melt = max(forcing["temperature_c"], 0) * 0.005 / 86400 * model.cos_beta
            snow_new = snow + snowfall_rate * dt
            melt_depth = min(snow_new, potential_melt * dt)
            snow_new -= melt_depth
            pet = forcing["pet_mm_h"] / 1000 / 3600 * model.cos_beta
            intercepted = min(rain_rate * dt, max(canopy_capacity - canopy, 0.0))
            canopy_new = canopy + intercepted
            canopy_evap = min(canopy_new, pet * dt)
            canopy_new -= canopy_evap
            root_pet = max(pet - canopy_evap / dt, 0.0)
            if scenario == "postfire":
                root_pet *= 0.35
            incoming = rain_rate - intercepted / dt + melt_depth / dt
            try:
                new_psi, new_surface, new_stage, info = model.advance(
                    psi, surface, stage, dt, incoming, root_pet
                )
            except RuntimeError:
                rejected += 1
                next_step = dt / 2.0
                if next_step < c.min_step_s:
                    raise RuntimeError(f"Adaptive Richards solve failed at hour {elapsed / 3600:.4f}")
                continue
            psi, surface, stage = new_psi, new_surface, new_stage
            snow, canopy = snow_new, canopy_new
            elapsed += dt
            total_substeps += 1
            min_accepted_step = min(min_accepted_step, dt)
            max_iterations = max(max_iterations, info["iterations"])
            max_local_residual = max(max_local_residual, info["local_max_residual_m3_s"])
            max_system_residual = max(max_system_residual, info["max_residual_m3_s"])
            max_accepted_pressure_increment = max(
                max_accepted_pressure_increment, info["final_pressure_increment_m"]
            )
            precipitation_volume = precip_rate * slope_area * dt
            et_volume = canopy_evap * slope_area + info["uptake_m3_s"] * dt
            outlet_volume = c.channel_area_m2 * (stage - model.channel_bed) / c.channel_tau_s * dt
            increments = {
                "precip_m3": precipitation_volume,
                "rain_m3": rain_rate * slope_area * dt,
                "melt_m3": melt_depth * slope_area,
                "et_m3": et_volume, "outlet_m3": outlet_volume,
                "runoff_m3": info["runoff_m3_s"] * dt,
                "gw_channel_m3": info["qtoe_m3_s"] * dt,
                "infiltration_m3": float(info["qtop_m3_s"].sum()) * dt,
                "substeps": 1, "iterations": info["iterations"],
            }
            for key, value in increments.items():
                hourly[key] += value
                if key in accumulated:
                    accumulated[key] += value
            next_step = min(c.max_step_s, dt * (1.35 if info["iterations"] < 20 else 1.05))
        save(hour + 1, hourly)
    final_stores = stores()
    total_change = sum(final_stores.values()) - initial_total
    residual = total_change - accumulated["precip_m3"] + accumulated["et_m3"] + accumulated["outlet_m3"]
    denominator = max(
        accumulated["precip_m3"], accumulated["et_m3"] + accumulated["outlet_m3"], 1.0e-12
    )
    runtime = time.perf_counter() - started
    budget = {
        "initial_stores_m3": initial_stores, "final_stores_m3": final_stores,
        **accumulated, "storage_change_m3": total_change,
        "signed_residual_m3": residual,
        "signed_residual_mm": residual / model.horizontal_area * 1000,
        "normalized_error": abs(residual) / denominator,
        "normalization_denominator_m3": denominator,
        "normalization_definition": "max(total precipitation, total ET + outlet volume); initial storage excluded",
        "max_abs_hourly_cumulative_residual_m3": max_hour_residual,
        "acceptance_tolerance_mm": 0.01,
        "accepted": abs(residual) / model.horizontal_area * 1000 < 0.01,
        "equation": "final_storage - initial_storage - precipitation + ET + outlet = residual",
    }
    return {
        "meta": {
            "schema_version": 1, "scenario": scenario,
            "model": "Synthetic 2D mixed-Richards hillslope reference; not coupled DHSVM/MODFLOW",
            "config": asdict(c), "units": "SI; displayed precipitation/SWE in mm; time in hours and UTC",
            "horizontal_crs": "none: synthetic local orthogonal slope coordinates",
            "vertical_datum": "arbitrary synthetic local elevation datum, metres upward",
            "coordinate_definition": "s downslope and d into plane; z=z0-s*sin(beta)-d*cos(beta); d is not vertical depth",
            "horizontal_area_m2": model.horizontal_area,
            "slope_area_m2": slope_area,
            "forcing_definition": "Synthetic hourly vertical precipitation/PET per horizontal projected area, converted by cos(beta) to slope-area flux",
            "equations": [
                "h = z + psi; Darcy q = -K(psi) grad(h); u = rho*g*psi",
                "d[theta(psi)+Ss*max(psi,0)]/dt = div(K grad(h)) - root_uptake",
                "Se=(1+(alpha*abs(psi))^n)^(-m), m=1-1/n for psi<0; Se=1 for psi>=0",
                "theta=theta_r+(theta_s-theta_r)*Se; K=Ks*sqrt(Se)*(1-(1-Se^(1/m))^m)^2",
                "dD/dt = throughfall + snowmelt - infiltration - max(D-D0,0)/tau_surface",
                "dV_channel/dt = bank_exchange + surface_runoff - V_channel/tau_channel",
            ],
            "boundaries": [
                "Impermeable divide and sloping base; no water crosses out of the section",
                "Top dry supply-limited flux / ponded pressure complementarity; implicit surface detention",
                "Toe bank above channel bed only: atmospheric seepage above stage and bidirectional exchange on submerged face fractions; toe below channel bed is no-flow; implicit reservoir",
            ],
            "assumptions": [
                "All geometry, materials, forcing and vegetation parameters are hypothetical; no field observations or assimilation",
                "Unified isotropic Richards continuum owns soil, vadose and groundwater storage; no duplicated shallow aquifer",
                "Temperature-index snowmelt (5 mm/deg C/day) conserves SWE; snow energy balance, glacier melt and frozen soil are omitted",
                "Uniform canopy bucket and Feddes-like uptake, oxygen-limited near saturation, use a single allocated PET demand",
                "Bare-soil evaporation, uptake compensation between depths and plant hydraulic redistribution are omitted; root_stress is an uptake multiplier (1 unstressed, 0 shut down)",
                "No retention hysteresis, preferential/macropore/fracture flow, density/heat transport or deep hydrothermal circulation",
                "Linear local surface detention and a lumped channel reservoir; no spatial overland routing or flood hydraulics",
                "The channel connects through the above-bed toe bank only (including partially wetted faces); horizontal channel-bed infiltration and regional groundwater flow beyond the closed toe wall are omitted",
                "Pressure diagnostic is antecedent hydrology only; no landslide, liquefaction, lateral spreading, fire or runout prediction",
                "Water-table depth is the deepest base-connected zero-pressure crossing; isolated perched saturation is not labelled regional water table",
                "Initial soil head is hydrostatic and channel starts with stored water; subsequent reservoir recession is an unspun-up synthetic transient",
                "Channel remains wet in the supported scenarios; complete river drying needs a supply-limited dry-bank boundary not implemented here, and causes an explicit failure",
                "The postfire option assumes Ks*0.3 in top 0.5 m, interception 0.2 mm and uptake demand*0.35; changes are uncalibrated and not universal fire effects",
            ],
            "parameter_provenance": "Illustrative synthetic values, not fitted or measured Washington soil parameters",
            "constitutive_parameters": {
                "soil": {"theta_r": 0.06, "theta_s": 0.43, "alpha_per_m": 1.2, "n": 1.6, "Ks_m_s": 1.5e-5},
                "weathered_substrate": {"theta_r": 0.06, "theta_s": 0.38, "alpha_per_m": 0.8, "n": 1.6, "Ks_m_s": 5e-6},
                "low_K_lens": {"theta_r": 0.06, "theta_s": 0.46, "alpha_per_m": 0.45, "n": 1.4, "Ks_m_s": 6e-8},
            },
            "sensor": {"x_m": float(model.x[model.sensor_ix]), "normal_depth_m": float(model.depth[model.sensor_iz]), "quantity": "rho*g*psi; negative values are matric suction"},
            "time_definition": "Artificial winter start 2026-01-01T00:00:00Z; end-of-hour states and preceding-hour mean rates",
        },
        "grid": {
            "x_m": model.x.tolist(), "depth_m": model.depth.tolist(),
            "elevation_m": model.z.tolist(), "surface_elevation_m": model.surface_z.tolist(),
            "layer_name": model.layer.tolist(), "dx_m": model.dx, "dz_m": model.dz,
            "width_m": c.width_m, "slope_rise_run": c.slope_rise_run,
        },
        "records": records, "snapshots": snapshots, "budget": budget,
        "diagnostics": {
            "wall_time_s": runtime, "simulated_hours": c.hours,
            "mean_wall_s_per_simulated_hour": runtime / c.hours,
            "accepted_substeps": total_substeps, "rejected_substeps": rejected,
            "min_accepted_step_s": min_accepted_step, "max_allowed_step_s": c.max_step_s,
            "max_nonlinear_iterations": max_iterations,
            "max_local_mass_residual_m3_s": max_local_residual,
            "max_system_residual_m3_s": max_system_residual,
            "local_residual_tolerance_m3_s": c.residual_tolerance_m3_s,
            "max_accepted_pressure_increment_m": max_accepted_pressure_increment,
            "pressure_increment_tolerance_m": c.pressure_tolerance_m,
            "spatial_convergence_tested": False,
            "forecast_or_ensemble_timing": False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("runs/toy_watershed/results.json"))
    parser.add_argument("--scenario", choices=["baseline", "intense_storm", "postfire"], default="baseline")
    parser.add_argument("--hours", type=int, default=120)
    parser.add_argument("--max-step-s", type=float, default=600)
    args = parser.parse_args(argv)
    result = run_toy(args.scenario, hours=args.hours, max_step_s=args.max_step_s)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({"output": str(args.output), "budget": result["budget"], "diagnostics": result["diagnostics"]}, indent=2))
    return 0 if result["budget"]["accepted"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
