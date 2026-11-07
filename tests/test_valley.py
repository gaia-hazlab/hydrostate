"""Synthetic 3D valley workflow checks, not field or hazard validation."""

import json
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

try:
    import numpy as np
    from hydrostate.valley import build_valley, run_valley, valley_boundaries
except ModuleNotFoundError as exc:
    if exc.name not in {"numpy", "scipy"}:
        raise
    np = None


@unittest.skipIf(np is None, "requires optional numpy/scipy numerics dependencies")
class ValleyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run_valley(hours=6)

    def test_real_cartesian_stairs_and_three_dimensional_exchange(self):
        result = self.result
        surface = np.array(result["grid"]["surface_elevation_m"])
        self.assertGreater(np.ptp(surface, axis=0).max(), 0)
        self.assertGreater(np.ptp(surface, axis=1).max(), 0)
        self.assertGreater(result["diagnostics"]["max_abs_internal_x_flux_m3_s"], 0)
        self.assertGreater(result["diagnostics"]["max_abs_internal_y_flux_m3_s"], 0)
        self.assertFalse(result["diagnostics"]["mesh_convergence_established"])

    def test_ledger_local_residual_and_nonponding_audit(self):
        result = self.result
        self.assertTrue(result["budget"]["accepted"])
        self.assertLess(abs(result["budget"]["signed_residual_mm"]), 1e-6)
        self.assertLess(result["diagnostics"]["max_local_residual_m3_s"], 1e-10)
        self.assertLess(result["diagnostics"]["max_implied_top_pressure_head_m"], 0)
        self.assertEqual(len(result["records"]), 7)

    def test_outlet_excludes_internal_east_facing_staircase_risers(self):
        model, psi, _, _ = build_valley()
        g = model.grid
        exposed = g.boundary_cells(2, 1)
        internal_risers = exposed[np.unravel_index(exposed, g.shape)[2] < g.shape[2] - 1]
        self.assertGreater(len(internal_risers), 0)
        regional_head, _ = valley_boundaries(model, .05)
        self.assertTrue(np.all(np.unravel_index(regional_head.cells, g.shape)[2] == g.shape[2] - 1))
        flux = model.face_fluxes(psi, boundaries=(regional_head,))
        self.assertFalse(np.isin(flux.boundary_cells, internal_risers).any())

    def test_json_null_mask_and_head_definition(self):
        result = self.result
        json.dumps(result, allow_nan=False)
        active = np.array(result["grid"]["active"])
        snapshot = result["snapshots"][-1]
        psi = np.array(snapshot["pressure_head_m"], dtype=float)
        head = np.array(snapshot["total_head_m"], dtype=float)
        z = np.array(result["grid"]["z_m"])[:, None, None]
        np.testing.assert_allclose((head - psi)[active], np.broadcast_to(z, active.shape)[active])
        self.assertTrue(np.isnan(psi[~active]).all())
        self.assertEqual(result["meta"]["outlet"]["id"], "east_valley_boundary")

    def test_identical_grid_refinement_comparison_keeps_all_active_cells(self):
        script = Path(__file__).resolve().parents[1] / "scripts" / "check_valley_refinement.py"
        spec = importlib.util.spec_from_file_location("valley_refinement_test", script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        comparison = module.compare(self.result, self.result)
        self.assertEqual(comparison["common_active_pressure_cells"], self.result["diagnostics"]["active_cells"])
        self.assertEqual(comparison["pressure_head_rmse_on_common_reference_centres_m"], 0.0)

    def test_ponding_incompatible_with_prescribed_infiltration_fails(self):
        with patch("hydrostate.valley.rainfall_mm_h", return_value=100.0):
            with self.assertRaisesRegex(ValueError, "requires ponding"):
                run_valley(hours=1)


if __name__ == "__main__":
    unittest.main()
