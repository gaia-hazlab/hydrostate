"""Targeted conservation/geometry/retry checks; analytical MMS lives separately."""

import unittest
from unittest.mock import patch

try:
    import numpy as np
    from hydrostate.flow import BoundaryFaces, ConvergenceError, Grid3D, RichardsModel, van_genuchten
except ModuleNotFoundError as exc:
    if exc.name not in {"numpy", "scipy"}:
        raise
    np = None


@unittest.skipIf(np is None, "requires optional numpy/scipy numerics dependencies")
class FlowTests(unittest.TestCase):
    def grid(self, active=None):
        return Grid3D((np.arange(4) + .5) * 2, (np.arange(3) + .5) * 1.5,
                      (np.arange(5) + .5) * .3, 2, 1.5, .3, active)

    def test_affine_total_head_and_darcy_flux_in_three_directions(self):
        g = self.grid()
        z, y, x = np.meshgrid(g.z_m, g.y_m, g.x_m, indexing="ij")
        affine = lambda xyz: 3 + .2 * xyz[:, 0] - .1 * xyz[:, 1] + .05 * xyz[:, 2]
        head = 3 + .2 * x - .1 * y + .05 * z
        psi = head - z
        model = RichardsModel(g, ks_m_s=2e-4, linear_capacity=.2)
        bc = []
        for axis in range(3):
            for side in (-1, 1):
                cells = g.boundary_cells(axis, side)
                bc.append(BoundaryFaces(cells, axis, side, "head", affine(g.face_centers_m(cells, axis, side))))
        residual = model.residual(psi, psi, 900, boundaries=bc)
        self.assertLess(np.max(abs(residual)), 1e-15)
        result = model.advance(psi, 900, boundaries=bc)
        np.testing.assert_allclose(result.pressure_head_m, psi, atol=1e-12)
        flux = result.fluxes
        differences = head.ravel()[flux.internal_i] - head.ravel()[flux.internal_j]
        # Unit gradients differ on all three axes, so exact affine flux tests
        # physical face areas/distances and z sign, not only zero equilibrium.
        for delta, spacing, area in ((1, g.dx_m, g.face_area_m2(2)),
                                     (4, g.dy_m, g.face_area_m2(1)),
                                     (12, g.dz_m, g.face_area_m2(0))):
            select = flux.internal_j - flux.internal_i == delta
            np.testing.assert_allclose(flux.internal_m3_s[select], 2e-4 * area / spacing * differences[select], atol=1e-16)

    def test_masked_heterogeneous_unsaturated_hydrostatic_equilibrium(self):
        active = np.ones((5, 3, 4), dtype=bool)
        active[3:, 1:, :2] = False
        g = self.grid(active)
        k = np.geomspace(1e-8, 1e-4, np.prod(g.shape)).reshape(g.shape)
        model = RichardsModel(g, ks_m_s=k)
        psi = np.where(g.active, 1.0 - g.elevation_m, np.nan)
        result = model.advance(psi, 3600)
        np.testing.assert_allclose(result.pressure_head_m, psi, atol=1e-12)
        self.assertLess(abs(result.water_budget_residual_m3), 1e-12)
        self.assertTrue(np.isnan(result.pressure_head_m[~active]).all())

    def test_internal_flux_cancellation_and_neumann_sign(self):
        g = self.grid()
        model = RichardsModel(g)
        rng = np.random.default_rng(4)
        psi = rng.uniform(-3, .5, g.shape)
        sources = rng.uniform(-1e-6, 1e-6, g.shape)
        west = g.boundary_cells(2, -1)
        bc = [BoundaryFaces(west, 2, -1, "flux", -2e-6)]
        residual = model.residual(psi, psi, 20, sources_m3_s=sources, boundaries=bc)
        expected = -2e-6 * len(west) - sources.sum()
        self.assertAlmostEqual(float(np.sum(residual)), expected, places=14)
        flux = model.face_fluxes(psi, boundaries=bc)
        self.assertTrue(np.all(flux.boundary_m3_s < 0))

    def test_failed_iteration_does_not_mutate_beginning_state(self):
        g = self.grid()
        model = RichardsModel(g, linear_capacity=.2)
        old = 1 - g.elevation_m
        before = old.copy()
        with self.assertRaises(ConvergenceError):
            model.advance(old, 100, sources_m3_s=1e-4, max_iterations=1)
        np.testing.assert_array_equal(old, before)
        result = model.advance(old, 100, sources_m3_s=1e-4)
        self.assertLess(abs(result.water_budget_residual_m3), 1e-10)

    def test_small_water_residual_does_not_bypass_pressure_correction(self):
        g = self.grid()
        model = RichardsModel(g, linear_capacity=1e-4)
        old = 1 - g.elevation_m
        result = model.advance(old, 100, sources_m3_s=1e-11 * g.cell_volume_m3)
        # Initial local residual is below the absolute water tolerance, but
        # the exact source-induced head change exceeds the pressure tolerance.
        np.testing.assert_allclose(result.pressure_head_m, old + 1e-5, atol=1e-12)
        self.assertLess(result.max_pressure_increment_m, 1e-7)

    def test_adaptive_rollback_preserves_exact_uniform_source_solution(self):
        g = self.grid()
        model = RichardsModel(g, linear_capacity=.2)
        old = 1 - g.elevation_m
        source = 1e-5 * g.cell_volume_m3
        original_advance = model.advance
        def reject_large(state, dt, **kwargs):
            if dt > 250:
                raise ConvergenceError("Injected failed solve to exercise rollback")
            return original_advance(state, dt, **kwargs)
        with patch.object(model, "advance", side_effect=reject_large):
            result = model.advance_interval(old, 1000, sources_m3_s=source, max_step_s=1000)
        self.assertGreaterEqual(result.rejected_steps, 2)
        np.testing.assert_allclose(result.pressure_head_m, old + 1e-5 * 1000 / .2, atol=1e-10)
        self.assertAlmostEqual(result.source_volume_m3, source * np.prod(g.shape) * 1000, places=11)
        self.assertLess(abs(result.water_budget_residual_m3), 1e-10)

    def test_nonlinear_closed_redistribution_bounds_and_mass(self):
        g = self.grid()
        model = RichardsModel(g)
        psi = np.full(g.shape, -1.0)
        result = model.advance_interval(psi, 1800, max_step_s=300)
        theta = model.theta(result.pressure_head_m)
        self.assertTrue(np.all((theta >= .06) & (theta <= .43)))
        self.assertLess(abs(result.water_budget_residual_m3), 1e-7)
        self.assertLess(result.pressure_head_m[-1].mean(), -1.0)
        self.assertGreater(result.pressure_head_m[0].mean(), -1.0)

    def test_boundary_and_grid_validation(self):
        g = self.grid()
        model = RichardsModel(g)
        cells = g.boundary_cells(2, -1)
        face = BoundaryFaces(cells, 2, -1, "flux", 0)
        with self.assertRaises(ValueError):
            model.face_fluxes(1 - g.elevation_m, boundaries=[face, face])
        with self.assertRaises(ValueError):
            Grid3D(np.array([0, 2, 5]), np.array([0]), np.array([0]), 2, 1, 1)
        with self.assertRaises(ValueError):
            model.advance(1 - g.elevation_m, 0)

    def test_van_genuchten_analytic_capacity_and_bounds(self):
        pressure = np.array([-100, -10, -1, -.1, 0, 1.0])
        theta, capacity, kr = van_genuchten(pressure)
        self.assertTrue(np.all(np.diff(theta) >= 0))
        self.assertTrue(np.all((theta >= .06) & (theta <= .43)))
        self.assertTrue(np.all((kr >= 0) & (kr <= 1)))
        epsilon = 1e-6
        finite_difference = (van_genuchten(pressure[:4] + epsilon)[0] - van_genuchten(pressure[:4] - epsilon)[0]) / (2 * epsilon)
        np.testing.assert_allclose(capacity[:4], finite_difference, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
