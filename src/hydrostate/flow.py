"""Three-dimensional orthogonal finite-volume mixed Richards reference.

Coordinates are Cartesian metres, z positive upward, and arrays are [z,y,x].
One shared volume flux per internal face enforces antisymmetric exchange.
For hydraulic head H=z+psi and outward volume flux Q, each active cell solves

    V [W(psi_new)-W(psi_old)] / dt + sum(Q_out) - source = 0.

For Richards flow W=theta(psi)+Ss*max(psi,0), using van Genuchten--Mualem.
The optional constant-K, W=C*psi mode is an explicitly linear mathematical
verification problem; its storage is not a physical soil water content.
Backward Euler is solved with conservative modified Picard iterations.
No nonlinear storage or conductivity is frozen between accepted steps.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable
import warnings

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve


@dataclass(frozen=True)
class Grid3D:
    """Uniform Cartesian cell centres and an immutable staircase active mask."""

    x_m: np.ndarray
    y_m: np.ndarray
    z_m: np.ndarray
    dx_m: float
    dy_m: float
    dz_m: float
    active: np.ndarray | None = None

    def __post_init__(self) -> None:
        for name, spacing in (("x_m", self.dx_m), ("y_m", self.dy_m), ("z_m", self.dz_m)):
            coordinate = np.array(getattr(self, name), dtype=float, copy=True)
            if (coordinate.ndim != 1 or coordinate.size < 1 or not np.all(np.isfinite(coordinate))
                    or not np.isfinite(spacing) or spacing <= 0):
                raise ValueError("Finite 1D centres and positive spacings are required")
            if coordinate.size > 1 and not np.allclose(np.diff(coordinate), spacing, rtol=1e-10, atol=1e-12):
                raise ValueError(f"{name} must be ascending and uniformly separated by its spacing")
            coordinate.setflags(write=False)
            object.__setattr__(self, name, coordinate)
        mask = np.ones(self.shape, dtype=bool) if self.active is None else np.array(self.active, dtype=bool, copy=True)
        if mask.shape != self.shape or not mask.any():
            raise ValueError("active must match (nz,ny,nx) and contain an active cell")
        mask.setflags(write=False)
        object.__setattr__(self, "active", mask)

    @property
    def shape(self) -> tuple[int, int, int]:
        return (len(self.z_m), len(self.y_m), len(self.x_m))

    @property
    def spacing_m(self) -> tuple[float, float, float]:
        return self.dz_m, self.dy_m, self.dx_m

    @property
    def elevation_m(self) -> np.ndarray:
        return np.broadcast_to(self.z_m[:, None, None], self.shape)

    @property
    def cell_volume_m3(self) -> float:
        return self.dx_m * self.dy_m * self.dz_m

    def face_area_m2(self, axis: int) -> float:
        if axis not in (0, 1, 2):
            raise ValueError("axis must be 0=z, 1=y or 2=x")
        return self.cell_volume_m3 / self.spacing_m[axis]

    def boundary_cells(self, axis: int, side: int) -> np.ndarray:
        """Full-array flat IDs of exposed faces, including staircase surfaces."""
        if axis not in (0, 1, 2) or side not in (-1, 1):
            raise ValueError("Use axis 0=z/1=y/2=x and side -1/+1")
        neighbor = np.roll(self.active, -side, axis=axis).copy()
        edge = [slice(None)] * 3
        edge[axis] = -1 if side == 1 else 0
        neighbor[tuple(edge)] = False
        return np.flatnonzero(self.active & ~neighbor)

    def face_centers_m(self, cells: np.ndarray, axis: int, side: int) -> np.ndarray:
        """Return xyz positions of selected face centres (one row per face)."""
        indices = np.unravel_index(np.asarray(cells, dtype=int), self.shape)
        xyz = np.column_stack((self.x_m[indices[2]], self.y_m[indices[1]], self.z_m[indices[0]]))
        xyz[:, 2 - axis] += side * self.spacing_m[axis] / 2
        return xyz


@dataclass(frozen=True)
class BoundaryFaces:
    """An exposed face group; omitted faces are impermeable.

    cells: full-array flat cell IDs. axis: 0=z,1=y,2=x. side: -1/+1.
    kind='head': TOTAL hydraulic head in metres at the boundary face.
    kind='flux': prescribed OUTWARD volume flux in m3/s per face (not m/s).
    Values may be scalar or one value per face. A physical face appears once.
    """

    cells: np.ndarray
    axis: int
    side: int
    kind: str
    values: np.ndarray | float


@dataclass
class FaceFluxes:
    """Signed shared internal i->j fluxes and outward boundary fluxes, m3/s."""

    internal_i: np.ndarray
    internal_j: np.ndarray
    internal_m3_s: np.ndarray
    boundary_cells: np.ndarray
    boundary_m3_s: np.ndarray


@dataclass
class StepResult:
    pressure_head_m: np.ndarray
    iterations: int
    max_local_residual_m3_s: float
    max_pressure_increment_m: float
    water_budget_residual_m3: float
    source_volume_m3: float
    boundary_out_volume_m3: float
    storage_change_m3: float
    fluxes: FaceFluxes


@dataclass
class IntervalResult:
    pressure_head_m: np.ndarray
    accepted_steps: int
    rejected_steps: int
    nonlinear_iterations: int
    min_step_s: float
    max_local_residual_m3_s: float
    max_pressure_increment_m: float
    water_budget_residual_m3: float
    source_volume_m3: float
    boundary_out_volume_m3: float
    storage_change_m3: float


class ConvergenceError(RuntimeError):
    """No accepted state was returned; caller can retry from its unchanged state."""


def van_genuchten(
    psi: np.ndarray, theta_r: np.ndarray | float = .06,
    theta_s: np.ndarray | float = .43, alpha_per_m: np.ndarray | float = 1.2,
    n: np.ndarray | float = 1.6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """theta, dtheta/dpsi [1/m], Mualem relative K with l=1/2; psi<0 suction."""
    p = np.asarray(psi, dtype=float)
    n = np.asarray(n)
    m = 1 - 1 / n
    a = np.asarray(alpha_per_m) * np.maximum(-p, 0)
    an = a ** n
    se = (1 + an) ** (-m)
    theta = theta_r + (np.asarray(theta_s) - theta_r) * se
    capacity = (np.asarray(theta_s) - theta_r) * m * n * alpha_per_m * a ** (n - 1) * (1 + an) ** (-m - 1)
    kr = np.sqrt(se) * (1 - (an / (1 + an)) ** m) ** 2
    return theta, capacity, kr


class RichardsModel:
    """Scalar isotropic conductivity, full 3D lateral/vertical Richards flow.

    Full-grid inputs/outputs keep inactive cells as NaN. Only active cells are
    solved. Materials accept scalar or grid-shaped arrays. Internal K uses a
    harmonic face mean; Dirichlet faces use cell K over half a cell distance.
    Thin fronts, sharp contrasts and stairs require explicit mesh refinement.
    """

    def __init__(
        self, grid: Grid3D, *, theta_r: np.ndarray | float = .06,
        theta_s: np.ndarray | float = .43, alpha_per_m: np.ndarray | float = 1.2,
        n: np.ndarray | float = 1.6, ks_m_s: np.ndarray | float = 1e-5,
        ss_per_m: np.ndarray | float = 1e-5,
        linear_capacity: np.ndarray | float | None = None,
    ):
        self.grid = grid
        self.ids = np.flatnonzero(grid.active)
        self.count = self.ids.size
        self._lookup = np.full(np.prod(grid.shape), -1, dtype=int)
        self._lookup[self.ids] = np.arange(self.count)
        self._z = grid.elevation_m.ravel()[self.ids]
        self.theta_r = self._material(theta_r)
        self.theta_s = self._material(theta_s)
        self.alpha_per_m = self._material(alpha_per_m)
        self.n = self._material(n)
        self.ks_m_s = self._material(ks_m_s)
        self.ss_per_m = self._material(ss_per_m)
        self.linear_capacity = None if linear_capacity is None else self._material(linear_capacity)
        if np.any(self.ks_m_s <= 0) or np.any(self.ss_per_m < 0):
            raise ValueError("Require Ks > 0 and Ss >= 0 in active cells")
        if self.linear_capacity is not None:
            if np.any(self.linear_capacity <= 0):
                raise ValueError("linear_capacity must be strictly positive")
        elif (np.any(self.theta_s <= self.theta_r) or np.any(self.theta_r < 0)
              or np.any(self.theta_s > 1) or np.any(self.alpha_per_m <= 0) or np.any(self.n <= 1)):
            raise ValueError("Invalid van Genuchten material")
        full_ids = np.arange(np.prod(grid.shape)).reshape(grid.shape)
        ii, jj, geometry = [], [], []
        for axis in (0, 1, 2):
            lo, hi = [slice(None)] * 3, [slice(None)] * 3
            lo[axis], hi[axis] = slice(None, -1), slice(1, None)
            lo, hi = tuple(lo), tuple(hi)
            valid = grid.active[lo] & grid.active[hi]
            ii.append(self._lookup[full_ids[lo][valid]])
            jj.append(self._lookup[full_ids[hi][valid]])
            geometry.append(np.full(int(valid.sum()), grid.face_area_m2(axis) / grid.spacing_m[axis]))
        self.edge_i, self.edge_j = np.concatenate(ii), np.concatenate(jj)
        self._edge_geometry = np.concatenate(geometry)

    def _material(self, value: np.ndarray | float) -> np.ndarray:
        result = np.broadcast_to(np.asarray(value, dtype=float), self.grid.shape).ravel()[self.ids].copy()
        if not np.all(np.isfinite(result)):
            raise ValueError("Active material/source values must be finite")
        return result

    def _unpack(self, value: np.ndarray) -> np.ndarray:
        result = np.asarray(value, dtype=float)
        if result.shape != self.grid.shape:
            raise ValueError(f"Expected grid shape {self.grid.shape}")
        compact = result.ravel()[self.ids]
        if not np.all(np.isfinite(compact)):
            raise ValueError("Active pressures must be finite")
        return compact

    def _expand(self, value: np.ndarray) -> np.ndarray:
        full = np.full(self.grid.shape, np.nan)
        full.ravel()[self.ids] = value
        return full

    def _properties(self, p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        if self.linear_capacity is not None:
            storage = self.linear_capacity * p
            return storage, self.linear_capacity, self.ks_m_s, storage
        theta, cap, kr = van_genuchten(p, self.theta_r, self.theta_s, self.alpha_per_m, self.n)
        storage = theta + self.ss_per_m * np.maximum(p, 0)
        capacity = cap + self.ss_per_m * (p >= 0)
        return storage, capacity, self.ks_m_s * kr, theta

    def storage(self, psi: np.ndarray) -> np.ndarray:
        return self._expand(self._properties(self._unpack(psi))[0])

    def capacity(self, psi: np.ndarray) -> np.ndarray:
        return self._expand(self._properties(self._unpack(psi))[1])

    def conductivity(self, psi: np.ndarray) -> np.ndarray:
        return self._expand(self._properties(self._unpack(psi))[2])

    def theta(self, psi: np.ndarray) -> np.ndarray:
        return self._expand(self._properties(self._unpack(psi))[3])

    def _boundaries(self, boundaries: Iterable[BoundaryFaces]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        cells, values, geometry, head = [], [], [], []
        seen: set[tuple[int, int, int]] = set()
        for bc in boundaries:
            selected = np.asarray(bc.cells, dtype=int).reshape(-1)
            if bc.kind not in ("head", "flux"):
                raise ValueError("Boundary kind is 'head' or 'flux'")
            valid = set(self.grid.boundary_cells(bc.axis, bc.side).tolist())
            if any(int(cell) not in valid for cell in selected):
                raise ValueError("Boundary cells must identify exposed active faces")
            for cell in selected:
                key = (int(cell), bc.axis, bc.side)
                if key in seen:
                    raise ValueError("Duplicate physical boundary face")
                seen.add(key)
            prescribed = np.broadcast_to(np.asarray(bc.values, dtype=float), selected.shape)
            if not np.all(np.isfinite(prescribed)):
                raise ValueError("Boundary values must be finite")
            cells.extend(self._lookup[selected])
            values.extend(prescribed)
            geometry.extend(np.full(selected.size, self.grid.face_area_m2(bc.axis) / (self.grid.spacing_m[bc.axis] / 2)))
            head.extend([bc.kind == "head"] * selected.size)
        return np.asarray(cells, dtype=int), np.asarray(values), np.asarray(geometry), np.asarray(head, dtype=bool)

    def _assemble(
        self, p: np.ndarray, old_storage: np.ndarray, dt_s: float,
        source: np.ndarray, boundary: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ):
        storage, capacity, conductivity, _ = self._properties(p)
        h = p + self._z
        i, j = self.edge_i, self.edge_j
        conductance = self._edge_geometry * 2 * conductivity[i] * conductivity[j] / np.maximum(conductivity[i] + conductivity[j], 1e-100)
        internal_flux = conductance * (h[i] - h[j])
        residual = (storage - old_storage) * self.grid.cell_volume_m3 / dt_s - source
        np.add.at(residual, i, internal_flux)
        np.add.at(residual, j, -internal_flux)
        diagonal = np.maximum(capacity, 1e-14) * self.grid.cell_volume_m3 / dt_s
        np.add.at(diagonal, i, conductance)
        np.add.at(diagonal, j, conductance)
        bc_cells, prescribed, bc_geometry, is_head = boundary
        bc_conductance = bc_geometry * conductivity[bc_cells] * is_head
        boundary_flux = np.where(is_head, bc_conductance * (h[bc_cells] - prescribed), prescribed)
        np.add.at(residual, bc_cells, boundary_flux)
        np.add.at(diagonal, bc_cells, bc_conductance)
        matrix = coo_matrix((
            np.r_[-conductance, -conductance, diagonal],
            (np.r_[i, j, np.arange(self.count)], np.r_[j, i, np.arange(self.count)]),
        ), shape=(self.count, self.count)).tocsr()
        fluxes = FaceFluxes(
            self.ids[i].copy(), self.ids[j].copy(), internal_flux,
            self.ids[bc_cells].copy(), boundary_flux,
        )
        return matrix, residual, fluxes, storage

    def residual(
        self, psi: np.ndarray, psi_old: np.ndarray, dt_s: float, *,
        sources_m3_s: np.ndarray | float = 0.0, boundaries: Iterable[BoundaryFaces] = (),
    ) -> np.ndarray:
        if not np.isfinite(dt_s) or dt_s <= 0:
            raise ValueError("dt_s must be positive and finite")
        p = self._unpack(psi)
        previous = self._properties(self._unpack(psi_old))[0]
        _, residual, _, _ = self._assemble(p, previous, dt_s, self._material(sources_m3_s), self._boundaries(boundaries))
        return self._expand(residual)

    def face_fluxes(self, psi: np.ndarray, *, boundaries: Iterable[BoundaryFaces] = ()) -> FaceFluxes:
        p = self._unpack(psi)
        return self._assemble(p, self._properties(p)[0], 1.0, np.zeros(self.count), self._boundaries(boundaries))[2]

    def advance(
        self, psi_old: np.ndarray, dt_s: float, *, sources_m3_s: np.ndarray | float = 0.0,
        boundaries: Iterable[BoundaryFaces] = (), tolerance_m3_s: float = 1e-10,
        head_tolerance_m: float = 1e-7, max_iterations: int = 80,
    ) -> StepResult:
        """Return an accepted state or raise; input arrays are never mutated.

        Both maximum LOCAL volume residual and pressure increment must pass.
        Storage is re-evaluated nonlinearly in every residual. Damping applies
        to trial iterations only; accepted moisture/storage is never clipped.
        """
        if (not np.isfinite(dt_s) or dt_s <= 0 or tolerance_m3_s <= 0
                or head_tolerance_m <= 0 or max_iterations < 1):
            raise ValueError("Positive time, tolerances and max_iterations required")
        p_old = self._unpack(psi_old)
        old_storage = self._properties(p_old)[0]
        p = p_old.copy()
        source = self._material(sources_m3_s)
        boundary = self._boundaries(boundaries)
        for iteration in range(1, max_iterations + 1):
            matrix, residual, fluxes, storage = self._assemble(p, old_storage, dt_s, source, boundary)
            maximum = float(np.max(abs(residual)))
            with warnings.catch_warnings():
                warnings.simplefilter("error", MatrixRankWarning)
                try:
                    increment = spsolve(matrix, -residual)
                except (MatrixRankWarning, RuntimeError) as exc:
                    raise ConvergenceError("Singular Richards linearization") from exc
            if not np.all(np.isfinite(increment)):
                raise ConvergenceError("Nonfinite Richards increment")
            biggest = float(np.max(abs(increment)))
            # A small absolute water residual can still imply a large head
            # correction when capacity/conductance are small. Check the actual
            # pending correction, including on the initial iterate.
            if maximum <= tolerance_m3_s and biggest <= head_tolerance_m:
                change = float((storage - old_storage).sum() * self.grid.cell_volume_m3)
                incoming = float(source.sum() * dt_s)
                outgoing = float(fluxes.boundary_m3_s.sum() * dt_s)
                return StepResult(self._expand(p), iteration, maximum, biggest,
                                  change + outgoing - incoming, incoming, outgoing, change, fluxes)
            damping = 1.0 if self.linear_capacity is not None else min(1.0, 1.0 / max(biggest, 1e-100))
            p += damping * increment
        raise ConvergenceError(f"No convergence in {max_iterations} iterations; local residual {maximum:g} m3/s")

    def advance_interval(
        self, psi_old: np.ndarray, duration_s: float, *,
        sources_m3_s: np.ndarray | float | Callable[[float], np.ndarray] = 0.0,
        boundaries: Iterable[BoundaryFaces] | Callable[[float], Iterable[BoundaryFaces]] = (),
        max_step_s: float = 600.0, min_step_s: float = .1,
        tolerance_m3_s: float = 1e-10, head_tolerance_m: float = 1e-7,
        on_accept: Callable[[StepResult, float], None] | None = None,
    ) -> IntervalResult:
        """Rollback/retry adaptive steps; callbacks receive interval-relative end seconds.

        Forcing callbacks must be deterministic and side-effect free because a
        rejected candidate can be requested again. Each accepted substep uses
        its own end-time source/head; there is no hidden accumulated state.
        on_accept(result, end_second) audits accepted substeps only; it must
        not mutate result. An audit exception fails the interval explicitly.
        """
        if (duration_s <= 0 or not np.isfinite(duration_s) or min_step_s <= 0
                or max_step_s < min_step_s or not np.isfinite(max_step_s)):
            raise ValueError("Require finite duration>0 and 0<min_step<=max_step")
        state = np.array(psi_old, copy=True)
        self._unpack(state)
        original = float(np.nansum(self.storage(state)) * self.grid.cell_volume_m3)
        elapsed, step = 0.0, min(max_step_s, duration_s)
        accepted = rejected = iterations = 0
        incoming = outgoing = max_residual = max_increment = 0.0
        minimum = duration_s
        static_boundaries = None if callable(boundaries) else tuple(boundaries)
        while elapsed < duration_s - 1e-10:
            dt = min(step, duration_s - elapsed)
            source = sources_m3_s(elapsed + dt) if callable(sources_m3_s) else sources_m3_s
            bc = boundaries(elapsed + dt) if callable(boundaries) else static_boundaries
            try:
                result = self.advance(state, dt, sources_m3_s=source, boundaries=bc,
                                      tolerance_m3_s=tolerance_m3_s, head_tolerance_m=head_tolerance_m)
            except ConvergenceError:
                rejected += 1
                step = dt / 2
                if step < min_step_s:
                    raise ConvergenceError(f"Adaptive step failed at interval second {elapsed:g}")
                continue
            if on_accept is not None:
                on_accept(result, elapsed + dt)
            state = result.pressure_head_m
            accepted += 1
            iterations += result.iterations
            incoming += result.source_volume_m3
            outgoing += result.boundary_out_volume_m3
            max_residual = max(max_residual, result.max_local_residual_m3_s)
            max_increment = max(max_increment, result.max_pressure_increment_m)
            minimum = min(minimum, dt)
            elapsed += dt
            step = min(max_step_s, dt * (1.4 if result.iterations < 15 else 1.1))
        change = float(np.nansum(self.storage(state)) * self.grid.cell_volume_m3) - original
        return IntervalResult(state, accepted, rejected, iterations, minimum, max_residual,
                              max_increment, change + outgoing - incoming, incoming, outgoing, change)
