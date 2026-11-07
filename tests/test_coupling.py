"""Actual-kernel fixed-volume mass reconciliation and hydraulic restart checks."""

from copy import deepcopy
import hashlib
import json
import unittest

try:
    import numpy as np
    from hydrostate.flow import BoundaryFaces, Grid3D, RichardsModel
    from hydrostate.coupling import (
        CouplingError, build_checkpoint, model_fingerprint, pressure_from_storage,
        reconcile_material_change, restore_checkpoint,
    )
except ModuleNotFoundError as exc:
    if exc.name not in {"numpy", "scipy"}:
        raise
    np = None


REFERENCE = {"horizontal_crs": "LOCAL:synthetic", "vertical_datum": "synthetic local elevation",
             "length_unit": "m", "z_positive": "up"}


def resign(payload):
    payload.pop("content_sha256", None)
    payload["content_sha256"] = hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    return payload


@unittest.skipIf(np is None, "requires optional numpy/scipy numerics dependencies")
class MaterialCouplingTests(unittest.TestCase):
    def setUp(self):
        active = np.ones((3, 2, 2), dtype=bool)
        active[2, 1, 1] = False
        self.grid = Grid3D([.5, 1.5], [.5, 1.5], [.25, .75, 1.25], 1, 1, .5, active)
        self.model = RichardsModel(self.grid, ss_per_m=1e-4)
        self.psi = np.where(active, 1 - self.grid.elevation_m, np.nan)

    def test_storage_inverse_spans_dry_unsaturated_and_saturated_pressure(self):
        psi = np.array([[[-20, -3], [-1, -.1]], [[0, .3], [2, 10]], [[-.01, -.5], [1, np.nan]]])
        recovered = pressure_from_storage(self.model, self.model.storage(psi))
        np.testing.assert_allclose(recovered, psi, rtol=1e-10, atol=1e-10, equal_nan=True)
        theta = self.model.theta(psi)
        self.assertGreater(self.model.storage(psi)[1, 1, 1], theta[1, 1, 1])

    def test_saturated_incompressible_branch_requires_an_explicit_pressure(self):
        model = RichardsModel(self.grid, ss_per_m=0)
        psi = np.where(self.grid.active, 2.0, np.nan)
        water = model.storage(psi)
        with self.assertRaisesRegex(CouplingError, "cannot determine pressure"):
            pressure_from_storage(model, water)
        np.testing.assert_array_equal(pressure_from_storage(model, water, saturated_pressure_head_m=psi), psi)
        with self.assertRaisesRegex(CouplingError, "nonnegative"):
            pressure_from_storage(model, water, saturated_pressure_head_m=-psi)
        with self.assertRaisesRegex(CouplingError, "exceeds storage capacity"):
            pressure_from_storage(model, water + np.where(self.grid.active, .001, np.nan), saturated_pressure_head_m=psi)

    def test_unrepresentable_dry_storage_and_pressure_envelope_are_rejected(self):
        for value in (.06, .05):
            with self.subTest(value=value), self.assertRaisesRegex(CouplingError, "no finite"):
                pressure_from_storage(self.model, np.where(self.grid.active, value, np.nan))
        high = np.where(self.grid.active, .43 + 1e-4 * 101, np.nan)
        with self.assertRaisesRegex(CouplingError, "admissible envelope"):
            pressure_from_storage(self.model, high)
        np.testing.assert_allclose(pressure_from_storage(self.model, high, pressure_bounds_m=(-1000, 200))[self.grid.active], 101)

    def test_ks_only_update_preserves_pressure_and_storage_exactly(self):
        target = RichardsModel(self.grid, ss_per_m=1e-4, ks_m_s=3e-5)
        original = self.psi.copy()
        result = reconcile_material_change(self.model, self.psi, target)
        np.testing.assert_array_equal(result.pressure_head_m, self.psi)
        np.testing.assert_array_equal(self.psi, original)
        self.assertEqual(result.absolute_residual_m3, 0)
        self.assertEqual(result.max_abs_cell_residual_m3, 0)

    def test_retention_and_storage_changes_preserve_each_cells_water(self):
        target = RichardsModel(self.grid, theta_s=.45, theta_r=.04, alpha_per_m=.8, n=1.5, ss_per_m=2e-4)
        water = self.model.storage(self.psi)
        before_model = model_fingerprint(self.model)
        result = reconcile_material_change(self.model, self.psi, target, stored_water=water)
        np.testing.assert_allclose(result.storage_after, water, atol=2e-15, rtol=0, equal_nan=True)
        self.assertLess(result.absolute_residual_m3, 1e-14)
        self.assertLess(result.max_abs_cell_residual_m3, 1e-15)
        self.assertEqual(model_fingerprint(self.model), before_model)
        self.assertGreater(float(np.nanmax(abs(result.pressure_head_m - self.psi))), .1)
        # theta can change despite exact W conservation because elastic storage
        # belongs in the hydraulic ledger but not physical liquid saturation.
        self.assertGreater(float(np.nanmax(abs(result.delta_theta))), 1e-6)

    def test_failed_material_update_is_atomic_and_cannot_hide_large_pressure(self):
        target = RichardsModel(self.grid, theta_s=.3, ss_per_m=1e-6)
        original_psi, original_model = self.psi.copy(), model_fingerprint(self.model)
        with self.assertRaisesRegex(CouplingError, "admissible envelope"):
            reconcile_material_change(self.model, self.psi, target)
        np.testing.assert_array_equal(self.psi, original_psi)
        self.assertEqual(model_fingerprint(self.model), original_model)

    def test_explicit_water_detects_material_mutation_and_new_grid_rejected(self):
        water = self.model.storage(self.psi)
        self.model.theta_s[:] += .01
        with self.assertRaisesRegex(CouplingError, "possible material mutation"):
            reconcile_material_change(self.model, self.psi, RichardsModel(self.grid), stored_water=water)
        moved = Grid3D([.51, 1.51], [.5, 1.5], [.25, .75, 1.25], 1, 1, .5, self.grid.active)
        with self.assertRaisesRegex(CouplingError, "conservative remeshing"):
            reconcile_material_change(self.model, self.psi, RichardsModel(moved))

    def test_linear_verification_storage_has_an_exact_inverse_but_is_not_mixed_with_water(self):
        model = RichardsModel(self.grid, linear_capacity=.2)
        np.testing.assert_allclose(pressure_from_storage(model, model.storage(self.psi)), self.psi, equal_nan=True)
        updated = RichardsModel(self.grid, linear_capacity=.4)
        result = reconcile_material_change(model, self.psi, updated)
        np.testing.assert_allclose(result.pressure_head_m, self.psi / 2, equal_nan=True)
        with self.assertRaisesRegex(CouplingError, "linear-verification"):
            reconcile_material_change(model, self.psi, self.model)

    def test_repeated_material_damage_and_recovery_does_not_accumulate_water(self):
        psi, model = self.psi, self.model
        original = model.storage(psi).copy()
        for index in range(20):
            fraction = .5 * (1 - np.cos(index * np.pi / 10))
            target = RichardsModel(self.grid, theta_s=.43 + .01 * fraction,
                                    alpha_per_m=1.2 * np.exp(.1 * fraction), ss_per_m=1e-4)
            update = reconcile_material_change(model, psi, target, stored_water=original)
            psi, model = update.pressure_head_m, target
            self.assertLess(update.max_abs_cell_residual_m3, 1e-15)
        np.testing.assert_allclose(model.storage(psi), original, atol=2e-15, rtol=0, equal_nan=True)


@unittest.skipIf(np is None, "requires optional numpy/scipy numerics dependencies")
class HydraulicCheckpointTests(unittest.TestCase):
    def setUp(self):
        self.grid = Grid3D([.5, 1.5], [.5, 1.5], [.25, .75], 1, 1, .5)
        self.model = RichardsModel(self.grid, ks_m_s=3e-6)
        self.psi = .2 - self.grid.elevation_m
        self.payload = build_checkpoint(self.model, self.psi, valid_time_utc="2026-01-02T02:00:00+02:00",
                                        issued_time_utc="2026-01-01T00:00:00Z", reference=REFERENCE,
                                        metadata={"forcing_version": "synthetic-v1"})

    def test_detached_roundtrip_complete_material_state_and_utc(self):
        state = restore_checkpoint(self.payload, expected_reference=REFERENCE, expected_model=self.model,
                                    expected_grid=self.grid, expected_valid_time_utc="2026-01-02T00:00:00Z")
        self.assertEqual(state.valid_time_utc, "2026-01-02T00:00:00Z")
        np.testing.assert_array_equal(state.pressure_head_m, self.psi)
        self.assertEqual(model_fingerprint(state.model), model_fingerprint(self.model))
        self.payload["metadata"]["forcing_version"] = "mutated-later"
        self.assertEqual(state.metadata["forcing_version"], "synthetic-v1")
        self.model.ks_m_s[:] *= 2
        self.assertFalse(np.array_equal(state.model.ks_m_s, self.model.ks_m_s))

    def test_actual_flow_restart_equivalence_including_material_update(self):
        top = self.grid.boundary_cells(0, 1)
        boundary = (BoundaryFaces(top, 0, 1, "flux", -1e-9),)
        first = self.model.advance_interval(self.psi, 900, boundaries=boundary, max_step_s=150)
        changed = RichardsModel(self.grid, theta_s=.44, ks_m_s=4e-6)
        update = reconcile_material_change(self.model, first.pressure_head_m, changed)
        uninterrupted = changed.advance_interval(update.pressure_head_m, 900, boundaries=boundary, max_step_s=150)
        checkpoint = build_checkpoint(changed, update.pressure_head_m, valid_time_utc="2026-01-02T00:15:00Z", reference=REFERENCE)
        restored = restore_checkpoint(checkpoint)
        resumed = restored.model.advance_interval(restored.pressure_head_m, 900, boundaries=boundary, max_step_s=150)
        np.testing.assert_array_equal(resumed.pressure_head_m, uninterrupted.pressure_head_m)
        self.assertEqual(resumed.storage_change_m3, uninterrupted.storage_change_m3)
        self.assertEqual(resumed.boundary_out_volume_m3, uninterrupted.boundary_out_volume_m3)

    def test_live_material_mutation_grid_clock_and_datum_mismatch_rejected(self):
        self.model.theta_s[:] += .001
        with self.assertRaisesRegex(CouplingError, "material fingerprint mismatch"):
            restore_checkpoint(self.payload, expected_model=self.model)
        with self.assertRaisesRegex(CouplingError, "clock mismatch"):
            restore_checkpoint(self.payload, expected_valid_time_utc="2026-01-03T00:00:00Z")
        with self.assertRaisesRegex(CouplingError, "reference/datum mismatch"):
            restore_checkpoint(self.payload, expected_reference={**REFERENCE, "vertical_datum": "other"})
        shifted = Grid3D([.51, 1.51], [.5, 1.5], [.25, .75], 1, 1, .5)
        with self.assertRaisesRegex(CouplingError, "grid/active-mask mismatch"):
            restore_checkpoint(self.payload, expected_grid=shifted)

    def test_content_units_masks_and_constitutive_inconsistency_rejected(self):
        corrupt = deepcopy(self.payload)
        corrupt["state"]["stored_water"]["values"][0][0][0] += .001
        with self.assertRaisesRegex(CouplingError, "content hash mismatch"):
            restore_checkpoint(corrupt)
        for section, name, units in (("state", "stored_water", "kg"), ("materials", "ks_m_s", "m/day")):
            corrupt = deepcopy(self.payload)
            corrupt[section][name]["units"] = units
            with self.subTest(section=section), self.assertRaisesRegex(CouplingError, "units mismatch"):
                restore_checkpoint(resign(corrupt))
        corrupt = deepcopy(self.payload)
        corrupt["state"]["theta"]["values"][0][0][0] += .001
        with self.assertRaisesRegex(CouplingError, "constitutive law"):
            restore_checkpoint(resign(corrupt))
        corrupt = deepcopy(self.payload)
        corrupt["grid"]["active"][0][0][0] = 1
        with self.assertRaisesRegex(CouplingError, "JSON booleans"):
            restore_checkpoint(resign(corrupt))
        corrupt = deepcopy(self.payload)
        corrupt["grid"]["dx_m"] = True
        with self.assertRaisesRegex(CouplingError, "numeric SI"):
            restore_checkpoint(resign(corrupt))

    def test_naive_clock_and_wrong_unit_reference_rejected(self):
        with self.assertRaisesRegex(CouplingError, "timezone-aware"):
            build_checkpoint(self.model, self.psi, valid_time_utc="2026-01-02T00:00:00", reference=REFERENCE)
        with self.assertRaisesRegex(CouplingError, "length_unit"):
            build_checkpoint(self.model, self.psi, valid_time_utc="2026-01-02T00:00:00Z",
                             reference={**REFERENCE, "length_unit": "ft"})

    def test_model_clock_rejects_silent_submicrosecond_truncation(self):
        for stamp in ("2026-01-02T00:00:00.1234567Z", "2026-01-02T00:00:00,123456789Z"):
            with self.subTest(stamp=stamp), self.assertRaisesRegex(CouplingError, "microsecond"):
                build_checkpoint(self.model, self.psi, valid_time_utc=stamp, reference=REFERENCE)
        payload = build_checkpoint(self.model, self.psi,
            valid_time_utc="2026-01-02T00:00:00.123456Z", reference=REFERENCE)
        self.assertEqual(payload["valid_time_utc"], "2026-01-02T00:00:00.123456Z")
        # The datahub integration CI installs pandas; the lightweight numerical
        # environment need not depend on it to reject overprecise ISO strings.
        try:
            import pandas as pd
        except ImportError:
            return
        with self.assertRaisesRegex(CouplingError, "nonzero nanoseconds"):
            build_checkpoint(self.model, self.psi,
                valid_time_utc=pd.Timestamp("2026-01-02T00:00:00.123456789Z"), reference=REFERENCE)
        payload = build_checkpoint(self.model, self.psi,
            valid_time_utc=pd.Timestamp("2026-01-02T00:00:00.123456Z"), reference=REFERENCE)
        self.assertEqual(payload["valid_time_utc"], "2026-01-02T00:00:00.123456Z")


if __name__ == "__main__":
    unittest.main()
