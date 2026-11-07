"""Analytical code verification of the actual hydrostate.flow kernel.

These checks distinguish equation/assembly accuracy from mass closure and
field validation. Manufactured sources are analytical continuum derivatives,
not residuals evaluated by the numerical operator being tested.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import time
from typing import Callable

import numpy as np
import scipy

from .flow import BoundaryFaces, Grid3D, RichardsModel, van_genuchten


def make_grid(n: int, *, lengths=(1.0, 1.2, 0.8), active=None) -> Grid3D:
    """A Cartesian cube of n**3 cells with deliberately unequal axis lengths."""
    lx, ly, lz = lengths
    dx, dy, dz = lx / n, ly / n, lz / n
    return Grid3D((np.arange(n) + 0.5) * dx,
                  (np.arange(n) + 0.5) * dy,
                  (np.arange(n) + 0.5) * dz, dx, dy, dz, active=active)


def coordinates(grid: Grid3D):
    z, y, x = np.meshgrid(grid.z_m, grid.y_m, grid.x_m, indexing="ij")
    return x, y, z


def head_boundaries(grid: Grid3D, head: Callable, t: float = 0.0):
    """Sample exact total head at each actual exposed face center."""
    x, y, z = coordinates(grid)
    spacing = (grid.dz_m, grid.dy_m, grid.dx_m)
    bc = []
    for axis in range(3):
        for side in (-1, 1):
            cells = np.asarray(grid.boundary_cells(axis, side), dtype=int)
            xyz = [x.ravel()[cells].copy(), y.ravel()[cells].copy(), z.ravel()[cells].copy()]
            xyz[2 - axis] += side * spacing[axis] / 2
            values = np.broadcast_to(head(*xyz, t), cells.shape).copy()
            bc.append(BoundaryFaces(cells=cells, axis=axis, side=side,
                                    kind="head", values=values))
    return tuple(bc)


def norms(error, active=None):
    """Volume-normalized L2 and Linf for equal-volume Cartesian cells."""
    e = np.asarray(error)[active] if active is not None else np.asarray(error).ravel()
    return {"l2_m": float(np.sqrt(np.mean(e * e))), "linf_m": float(np.max(np.abs(e)))}


def convergence_rates(records, scale_key: str):
    rates = []
    for coarse, fine in zip(records, records[1:]):
        ratio = coarse[scale_key] / fine[scale_key]
        rates.append({
            "coarse": coarse[scale_key], "fine": fine[scale_key],
            "l2_order": float(np.log(coarse["l2_m"] / fine["l2_m"]) / np.log(ratio)),
            "linf_order": float(np.log(coarse["linf_m"] / fine["linf_m"]) / np.log(ratio)),
        })
    return rates


def analytic_vg(psi, *, theta_r=0.06, theta_s=0.43, alpha=1.2, n=1.6, ks=1e-5):
    """Independent analytical values and chain-rule derivatives, psi < 0 only.

    Source generation for nonlinear MMS uses these continuum formulas rather
    than hydrostate.flow's retention or discrete flux/residual routines.
    """
    p = np.asarray(psi, dtype=float)
    if np.any(p >= 0):
        raise ValueError("The smooth analytical VG oracle requires psi < 0")
    m = 1 - 1 / n
    a = -alpha * p
    se = (1 + a ** n) ** (-m)
    dse = m * n * alpha * a ** (n - 1) * (1 + a ** n) ** (-m - 1)
    b = 1 - (1 - se ** (1 / m)) ** m
    db_dse = (1 - se ** (1 / m)) ** (m - 1) * se ** (1 / m - 1)
    theta = theta_r + (theta_s - theta_r) * se
    capacity = (theta_s - theta_r) * dse
    k = ks * se ** 0.5 * b ** 2
    dk = ks * (0.5 * se ** (-0.5) * b ** 2
               + 2 * se ** 0.5 * b * db_dse) * dse
    return theta, capacity, k, dk


def constitutive_check():
    psi = np.array([-30.0, -10.0, -2.0, -0.5, -0.1])
    theta, capacity, kr = van_genuchten(psi)
    a_theta, a_capacity, k, dk = analytic_vg(psi)
    epsilon = 1e-6
    _, _, kp = van_genuchten(psi + epsilon)
    _, _, km = van_genuchten(psi - epsilon)
    dk_fd = (kp - km) * 1e-5 / (2 * epsilon)
    theta_error = float(np.max(abs(theta - a_theta)))
    capacity_error = float(np.max(abs(capacity - a_capacity)))
    k_relative_error = float(np.max(abs(kr * 1e-5 - k) / k))
    dk_relative_error = float(np.max(abs(dk_fd - dk) / abs(dk)))
    # Verify the *model's* storage and capacity, including its elastic branch.
    grid = make_grid(2)
    model = RichardsModel(grid, ss_per_m=2e-5)
    storage_errors = []
    for p in (-2.0, -0.1, 0.1, 2.0):
        state = np.full(grid.elevation_m.shape, p)
        derivative = (model.storage(state + epsilon) - model.storage(state - epsilon)) / (2 * epsilon)
        storage_errors.append(float(np.max(abs(derivative - model.capacity(state)))))
    passed = (theta_error < 1e-13 and capacity_error < 1e-13
              and k_relative_error < 1e-10 and dk_relative_error < 2e-5
              and max(storage_errors) < 1e-8)
    return {"passed": passed, "theta_max_absolute_error": theta_error,
            "capacity_max_absolute_error_per_m": capacity_error,
            "conductivity_max_relative_error": k_relative_error,
            "conductivity_derivative_max_relative_error": dk_relative_error,
            "storage_derivative_max_absolute_error_per_m": max(storage_errors),
            "derivative_scope": "Smooth unsaturated points and separate positive-head elastic branch; not the branch junction"}


def linear_head_check():
    grid = make_grid(5, lengths=(2.5, 1.5, 1.0))
    x, y, z = coordinates(grid)
    k, c, dt = 1e-3, 0.2, 20.0
    def exact(x, y, z, t):
        return 4.0 + 0.3 * x - 0.2 * y + 0.4 * z
    psi = exact(x, y, z, 0) - z
    model = RichardsModel(grid, ks_m_s=k, linear_capacity=c)
    bc = head_boundaries(grid, exact)
    residual = model.residual(psi, psi, dt, boundaries=bc)
    result = model.advance(psi, dt, boundaries=bc,
                           tolerance_m3_s=1e-12, head_tolerance_m=1e-10)
    error = norms(result.pressure_head_m - psi)
    # Test each externally imposed flux, including half-cell distances/signs.
    # Boundary flux is independently inferred from each face's physical Darcy law.
    gradients = (0.4, -0.2, 0.3)
    spacing = (grid.dz_m, grid.dy_m, grid.dx_m)
    cell_volume = float(grid.cell_volume_m3)
    expected_faces = []
    for b in bc:
        expected_faces.extend(np.full(len(b.cells), -k * gradients[b.axis] * b.side
                                      * cell_volume / spacing[b.axis]))
    fluxes = model.face_fluxes(psi, boundaries=bc)
    boundary_error = float(np.max(abs(fluxes.boundary_m3_s - expected_faces)))
    locations = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    displacement = locations[fluxes.internal_j] - locations[fluxes.internal_i]
    distance = np.linalg.norm(displacement, axis=1)
    expected_internal = (-k * (displacement @ np.array([0.3, -0.2, 0.4]))
                         * cell_volume / distance**2)
    internal_error = float(np.max(abs(fluxes.internal_m3_s - expected_internal)))
    return {"passed": (error["linf_m"] < 1e-10
                       and float(np.nanmax(abs(residual))) < 1e-12
                       and boundary_error < 1e-12 and internal_error < 1e-12),
            **error, "max_local_residual_m3_s": float(np.nanmax(abs(residual))),
            "max_boundary_flux_error_m3_s": boundary_error,
            "max_internal_flux_error_m3_s": internal_error,
            "analytic_boundary_net_out_m3_s": float(np.sum(expected_faces)),
            "analytic_darcy_flux_xyz_m_s": [-0.3 * k, 0.2 * k, -0.4 * k],
            "scope": "Affine 3-D total head, unequal spacings, exact half-cell head boundaries and constant conductivity"}


def prescribed_source_flux_check():
    grid = make_grid(1)
    model = RichardsModel(grid, ks_m_s=0.01, linear_capacity=0.2)
    old = np.full(grid.shape, 2.0)
    cell = grid.boundary_cells(2, 1)
    bc = (BoundaryFaces(cell, 2, 1, "flux", 0.006),)
    source, dt = 0.002, 3.0
    exact = old + (source - 0.006) * dt / (0.2 * grid.cell_volume_m3)
    result = model.advance(old, dt, sources_m3_s=source, boundaries=bc,
                           tolerance_m3_s=1e-12, head_tolerance_m=1e-10)
    error = float(np.max(abs(result.pressure_head_m - exact)))
    return {"passed": error < 1e-10 and abs(result.water_budget_residual_m3) < 1e-12,
            "maximum_pressure_error_m": error,
            "source_volume_m3": result.source_volume_m3,
            "outward_boundary_volume_m3": result.boundary_out_volume_m3,
            "storage_change_m3": result.storage_change_m3,
            "scope": "Single-cell exact linear storage under positive source and positive outward prescribed volume flux; tests units and signs independently of Darcy head boundaries"}


def hydrostatic_mask_check():
    n = 5
    active = np.ones((n, n, n), dtype=bool)
    active[-2:, 1:3, :2] = False
    grid = make_grid(n, active=active)
    _, _, z = coordinates(grid)
    rng = np.random.default_rng(817)
    ks = 10 ** rng.uniform(-7, -3, size=z.shape)
    model = RichardsModel(grid, ks_m_s=ks)
    psi = np.where(active, 0.43 - z, np.nan)
    residual = model.residual(psi, psi, 100.0)
    result = model.advance(psi, 100.0, tolerance_m3_s=1e-12, head_tolerance_m=1e-10)
    error = norms(result.pressure_head_m - psi, active)
    return {"passed": error["linf_m"] < 1e-10 and float(np.nanmax(abs(residual))) < 1e-12,
            **error, "max_local_residual_m3_s": float(np.nanmax(abs(residual))),
            "active_cells": int(active.sum()), "inactive_cells": int((~active).sum()),
            "conductivity_range_m_s": [float(ks.min()), float(ks.max())],
            "scope": "Nonlinear retention across saturation, heterogeneous K, stair-step mask, all exposed faces sealed"}


def no_flow_redistribution_check():
    grid = make_grid(4)
    x, y, z = coordinates(grid)
    model = RichardsModel(grid, ks_m_s=1e-4)
    psi = -1.0 + 0.07 * np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y / 1.2)
    result = model.advance(psi, 10.0, tolerance_m3_s=1e-12, head_tolerance_m=1e-9)
    change = float(np.sum(model.storage(result.pressure_head_m) - model.storage(psi))
                   * grid.cell_volume_m3)
    moved = float(np.max(abs(result.pressure_head_m - psi)))
    return {"passed": abs(change) < 1e-9 and moved > 1e-6,
            "signed_storage_change_m3": change, "max_pressure_change_m": moved,
            "max_local_residual_m3_s": result.max_local_residual_m3_s,
            "scope": "Non-equilibrium nonlinear 3-D closed-domain redistribution; conservation is necessary but not a front-accuracy benchmark"}


def linear_jacobian_check():
    """Residual directional derivatives must form an SPD linear diffusion system.

    This checks the actual residual independently of its solver's matrix. It
    does not incorrectly demand that a frozen-conductivity Picard matrix equal
    the full nonlinear Newton Jacobian.
    """
    grid = make_grid(2, lengths=(1.0, 1.4, 0.8))
    k, c, dt = 0.03, 0.2, 0.7
    model = RichardsModel(grid, ks_m_s=k, linear_capacity=c)
    _, _, z = coordinates(grid)
    psi = 3.0 - z
    base = model.residual(psi, psi, dt).ravel()
    columns = []
    for j in range(psi.size):
        delta = np.zeros_like(psi)
        delta.ravel()[j] = 1e-5
        columns.append((model.residual(psi + delta, psi, dt).ravel() - base) / 1e-5)
    jac = np.column_stack(columns)
    compact = model._unpack(psi)
    solver_matrix = model._assemble(compact, model._properties(compact)[0], dt,
                                    np.zeros(model.count), model._boundaries(()))[0].toarray()
    assembly_error = float(np.max(abs(jac - solver_matrix)))
    symmetry = float(np.max(abs(jac - jac.T)))
    eigenvalue = float(np.linalg.eigvalsh(0.5 * (jac + jac.T)).min())
    row_sum_error = float(np.max(abs(jac.sum(axis=1) - c * grid.cell_volume_m3 / dt)))
    off_diagonal = jac.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    offmax = float(off_diagonal.max())
    return {"passed": symmetry < 1e-9 and eigenvalue > 0 and row_sum_error < 1e-9
                      and offmax < 1e-10 and assembly_error < 1e-9,
            "symmetry_max_error_m2_s": symmetry, "minimum_eigenvalue_m2_s": eigenvalue,
            "storage_row_sum_max_error_m2_s": row_sum_error,
            "largest_off_diagonal_m2_s": offmax,
            "actual_solver_matrix_derivative_max_error_m2_s": assembly_error,
            "scope": "Linear residual Jacobian: positive storage, shared symmetric diffusion, nonpositive neighbor couplings; no-flow boundaries"}


def linear_spatial_mms(ns=(5, 9, 17)):
    """All-three-axis sinusoid, linear time amplitude: backward Euler is exact in time."""
    lx, ly, lz = lengths = (1.0, 1.2, 0.8)
    kx, ky, kz = np.pi / lx, np.pi / ly, np.pi / lz
    wavenumber_sq = kx * kx + ky * ky + kz * kz
    k, c, amplitude, final_time, time_scale_s = 0.05, 0.25, 0.2, 0.3, 1.0
    def shape(x, y, z):
        return np.sin(kx * x) * np.sin(ky * y) * np.sin(kz * z)
    def head(x, y, z, t):
        return 3.0 + amplitude * (1 + t / time_scale_s) * shape(x, y, z)
    records = []
    for n in ns:
        grid = make_grid(n, lengths=lengths)
        x, y, z = coordinates(grid)
        old = head(x, y, z, 0) - z
        source = amplitude * (c / time_scale_s + k * wavenumber_sq
                              * (1 + final_time / time_scale_s)) * shape(x, y, z)
        model = RichardsModel(grid, ks_m_s=k, linear_capacity=c)
        result = model.advance(old, final_time,
                               sources_m3_s=source * grid.cell_volume_m3,
                               boundaries=head_boundaries(grid, head, final_time),
                               tolerance_m3_s=1e-12, head_tolerance_m=1e-10)
        records.append({"n": n, "h_m": max(grid.dx_m, grid.dy_m, grid.dz_m),
                        **norms(result.pressure_head_m - (head(x, y, z, final_time) - z)),
                        "max_local_residual_m3_s": result.max_local_residual_m3_s,
                        "budget_residual_m3": result.water_budget_residual_m3})
    rates = convergence_rates(records, "h_m")
    return {"passed": all(1.8 < r["l2_order"] < 2.3 and 1.8 < r["linf_order"] < 2.3 for r in rates),
            "records": records, "rates": rates, "expected_order": 2,
            "exact_total_head": "3 m + (0.2 m)*(1+t/(1 s))*sin(pi*x/(1 m))*sin(pi*y/(1.2 m))*sin(pi*z/(0.8 m))",
            "continuum_source_per_s": "A*[C/tau + K*(kx²+ky²+kz²)*(1+t/tau)]*P(x,y,z); tau=1 s",
            "scope": "Constant-K, linear-storage branch of the SAME 3-D flow kernel; midpoint source, cell-center norms; odd grids align the exact field maximum across meshes, avoiding Linf sampling artifacts; no time-truncation error for this linear-in-time exact solution"}


def temporal_mms(steps=(4, 8, 16)):
    """Uniform-in-space total head removes spatial truncation from the time test."""
    grid = make_grid(5)
    x, y, z = coordinates(grid)
    c, k, final_time, amplitude, decay = 0.25, 0.02, 0.4, 0.3, 2.0
    def head(x, y, z, t):
        return np.zeros_like(x) + 3.0 + amplitude * np.exp(-decay * t)
    records = []
    for count in steps:
        model = RichardsModel(grid, ks_m_s=k, linear_capacity=c)
        psi = head(x, y, z, 0) - z
        dt = final_time / count
        max_residual = 0.0
        for step in range(1, count + 1):
            t = step * dt
            source = np.full(z.shape, -c * decay * amplitude * np.exp(-decay * t)
                             * grid.cell_volume_m3)
            result = model.advance(psi, dt, sources_m3_s=source,
                                   boundaries=head_boundaries(grid, head, t),
                                   tolerance_m3_s=1e-12, head_tolerance_m=1e-10)
            psi = result.pressure_head_m
            max_residual = max(max_residual, result.max_local_residual_m3_s)
        records.append({"steps": count, "dt_s": dt,
                        **norms(psi - (head(x, y, z, final_time) - z)),
                        "max_local_residual_m3_s": max_residual})
    rates = convergence_rates(records, "dt_s")
    return {"passed": all(0.85 < r["l2_order"] < 1.15 and 0.85 < r["linf_order"] < 1.15 for r in rates),
            "records": records, "rates": rates, "expected_order": 1,
            "exact_total_head": "3 + 0.3*exp(-2*t)",
            "continuum_source_per_s": "-C*2*0.3*exp(-2*t)",
            "scope": "Spatially uniform total head has exact FV Laplacian zero; measures backward-Euler temporal error with time-dependent exact head boundaries"}


def nonlinear_spatial_mms(ns=(5, 9, 17), *, dt_s=1.0):
    """Smooth unsaturated VG field with analytical K' gravity/gradient source.

    Stationarity removes time-discretization error. The forcing is not formed
    by applying the numerical operator to the target field.
    """
    lx, ly, lz = (1.0, 1.2, 0.8)
    kx, ky, kz = np.pi / lx, np.pi / ly, np.pi / lz
    amplitude, psi_base, conductivity, dt = 0.12, -0.8, 0.05, dt_s
    def exact_psi(x, y, z):
        return psi_base + amplitude * np.sin(kx*x) * np.sin(ky*y) * np.sin(kz*z)
    def head(x, y, z, t):
        return exact_psi(x, y, z) + z
    records = []
    for n in ns:
        grid = make_grid(n)
        x, y, z = coordinates(grid)
        psi = exact_psi(x, y, z)
        sx, sy, sz = np.sin(kx*x), np.sin(ky*y), np.sin(kz*z)
        px = amplitude * kx * np.cos(kx*x) * sy * sz
        py = amplitude * ky * sx * np.cos(ky*y) * sz
        pz = amplitude * kz * sx * sy * np.cos(kz*z)
        lap = -(kx*kx + ky*ky + kz*kz) * (psi - psi_base)
        _, _, kval, dk = analytic_vg(psi, ks=conductivity)
        source = -dk * (px*px + py*py + pz*pz + pz) - kval * lap
        model = RichardsModel(grid, ks_m_s=conductivity)
        result = model.advance(psi, dt, sources_m3_s=source * grid.cell_volume_m3,
                               boundaries=head_boundaries(grid, head),
                               tolerance_m3_s=1e-12, head_tolerance_m=1e-9,
                               max_iterations=120)
        records.append({"n": n, "h_m": max(grid.dx_m, grid.dy_m, grid.dz_m),
                        **norms(result.pressure_head_m - psi),
                        "max_local_residual_m3_s": result.max_local_residual_m3_s,
                        "budget_residual_m3": result.water_budget_residual_m3})
    rates = convergence_rates(records, "h_m")
    return {"passed": None, "acceptance_required": False,
            "status": "diagnostic; nonlinear asymptotic Linf order not established",
            "asymptotic_linf_verified": False,
            "records": records, "rates": rates,
            "trial_duration_s": dt_s,
            "acceptance": "Diagnostic only: report measured errors/orders; no formal nonlinear asymptotic-order acceptance claimed",
            "exact_pressure_head": "-0.8 + 0.12*sin(pi*x/1)*sin(pi*y/1.2)*sin(pi*z/0.8)",
            "continuum_source_per_s": "-K'(psi)*(|grad psi|² + dpsi/dz) - K(psi)*laplacian(psi)",
            "scope": "Stationary smooth unsaturated nonlinear VG balance including gravity; does not validate sharp fronts, saturation transitions, hysteresis or heterogeneous-interface accuracy"}


def run_verification():
    started = time.perf_counter()
    checks = {
        "constitutive_derivatives": constitutive_check(),
        "affine_darcy_head": linear_head_check(),
        "prescribed_source_and_flux": prescribed_source_flux_check(),
        "heterogeneous_masked_hydrostatic": hydrostatic_mask_check(),
        "nonlinear_closed_redistribution": no_flow_redistribution_check(),
        "linear_residual_jacobian": linear_jacobian_check(),
        "linear_3d_spatial_mms": linear_spatial_mms(),
        "backward_euler_temporal_mms": temporal_mms(),
        "nonlinear_3d_spatial_mms": nonlinear_spatial_mms(),
        "nonlinear_short_duration_diagnostic": nonlinear_spatial_mms(dt_s=0.1),
    }
    kernel = Path(__file__).with_name("flow.py")
    return {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
            "passed": all(check["passed"] for check in checks.values()
                          if check.get("acceptance_required", True)),
            "acceptance_scope": "Required checks pass analytical linear spatial/time order, local constitutive/assembly checks and hydrostatic/conservation invariants. Nonlinear MMS remains an explicitly unaccepted asymptotic-order diagnostic.",
            "kernel": "hydrostate.flow.RichardsModel (same advance/residual as valley runs)",
            "kernel_sha256": hashlib.sha256(kernel.read_bytes()).hexdigest(),
            "verification_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__, "scipy_version": scipy.__version__,
            "wall_time_s": time.perf_counter()-started,
            "checks": checks,
            "limitations": [
                "Code verification is not field validation or uncertainty calibration.",
                "No proof of correctness for arbitrary meshes, anisotropy, nonsmooth coefficients or boundary conditions.",
                "Nonlinear infiltration fronts, capillary barriers, moving water tables, surface coupling and external reference comparisons remain separate validation/verification tasks.",
                "No DHSVM/MODFLOW interoperability, DA skill, ecology/mechanics or operational warning claim follows from this suite.",
            ],
            "references": [
                {"title": "Celia et al. 1990: mass-conservative mixed Richards formulation", "url": "https://doi.org/10.1029/WR026i007p01483"},
                {"title": "van Genuchten 1980: retention and conductivity", "url": "https://doi.org/10.2136/sssaj1980.03615995004400050002x"},
                {"title": "Salari and Knupp 2000: method of manufactured solutions", "url": "https://doi.org/10.2172/759450"},
                {"title": "NASA: Verification of Software: The Textbook and Real Problems", "url": "https://ntrs.nasa.gov/citations/20080010138"},
            ]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("runs/verification/report.json"))
    args = parser.parse_args(argv)
    report = run_verification()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "passed": report["passed"],
                      "checks": {key: check["passed"] for key, check in report["checks"].items()},
                      "wall_time_s": report["wall_time_s"]}, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
