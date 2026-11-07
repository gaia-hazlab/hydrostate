"""Typed interoperability, distinct model times, and actual Gaia round trips."""
from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

try:
    import numpy as np
    import xarray as xr
    from hydrostate.datahub import (CONTRACT, EVENT_GROUP, GAIA_COMMIT, INTERVAL_GROUP,
        MATERIAL_GROUP, STATE_GROUP, STATIC_GROUP, Field, build_tree, open_tree,
        save_tree, source_attrs, validate_tree, checkpoint_from_tree, canonical_reference)
except ImportError:
    np = xr = None


@unittest.skipIf(xr is None, "Install Hydrostate datahub dependencies")
class DatahubTests(unittest.TestCase):
    def arguments(self, *, projected=False):
        active = np.ones((2, 2, 2), dtype=bool)
        active[1, 1, 1] = False
        offset = 500_000.0 if projected else 0.0
        grid = SimpleNamespace(x_m=offset + np.array([.5, 1.5]),
            y_m=(5_200_000 if projected else 0) + np.array([.5, 1.5]),
            z_m=np.array([.5, 1.5]), dx_m=1., dy_m=1., dz_m=1., active=active)
        theta = np.full((2, 2, 2, 2), .25, dtype=np.float32)
        theta[:, ~active] = np.nan
        psi = np.full(theta.shape, -.8, dtype=np.float64)
        psi[:, ~active] = np.nan
        facies = np.full(active.shape, 2, dtype=np.int16)
        facies[~active] = -1
        material = np.full(theta.shape, .4, dtype=np.float64)
        material[:, ~active] = np.nan
        return dict(grid=grid,
            reference={"horizontal_crs": "EPSG:32610" if projected else "LOCAL_SYNTHETIC_CARTESIAN",
                       "vertical_datum": "NAVD88" if projected else "synthetic floor, metres",
                       "length_unit": "m", "z_positive": "up"},
            provenance=source_attrs(title="Synthetic coupled test", description="No field validation",
                source_commit="test-sha", lineage="Constructed synthetic fields"),
            state_times_utc=["2026-01-01T00:00:00Z", "2026-01-01T01:00:00Z"],
            state_fields={"pressure_head_m": Field(psi, "m"), "water_content": Field(theta, "m3/m3"),
                          "total_water_storage_fraction": Field(theta.astype(np.float64), "m3/m3")},
            static_fields={"facies_id": Field(facies, "1", {"missing_value_code": -1}),
                "facies_label": Field(np.full(active.shape, "volcanic alluvium with long text"), "1")},
            material_fields={"theta_s": Field(material, "m3/m3")},
            material_epochs=["initial", "quake-1"],
            interval_bounds_utc=[["2026-01-01T00:00:00Z", "2026-01-01T01:00:00Z"]],
            interval_fields={"boundary_out_volume_m3": Field(np.array([0.], dtype=np.float64), "m3", {"cell_methods": "time: sum"}),
                             "source_rate_m3_s": Field(np.array([1e-6], dtype=np.float32), "m3/s", {"cell_methods": "time: mean"})},
            event_times_utc=["2026-01-01T00:00:37.123456789Z"], event_ids=["quake-1"],
            event_fields={"kind": Field(np.array(["earthquake-conditioned material jump"]), "1"),
                "material_epoch_before": Field(np.array(["initial"]), "1"),
                "material_epoch_after": Field(np.array(["quake-1"]), "1"),
                "water_budget_residual_m3": Field(np.array([0.]), "m3"),
                "pressure_head_before_m": Field(psi[:1], "m"),
                "pressure_head_after_m": Field(psi[:1] - .1, "m")},
            checkpoint_payload={"schema": "test.coupled/1", "valid_time_utc": "2026-01-01T01:00:00Z",
                "recovery_modes": [{"amplitude": .23, "tau_s": 86400.}, {"amplitude": .12, "tau_s": 8640000.}],
                "complete_memory_text": "memory α" * 10000})

    def test_canonical_dtypes_and_local_reference(self):
        tree = build_tree(**self.arguments())
        self.assertEqual(tree.attrs["hydrostate:contract"], CONTRACT)
        self.assertEqual(tree.attrs["hydrostate:gaia-contract-commit"], GAIA_COMMIT)
        state = tree[STATE_GROUP].to_dataset()
        self.assertEqual(state.water_content.dims, ("time", "z", "y", "x"))
        self.assertEqual(state.water_content.dtype, np.dtype("float32"))
        self.assertEqual(state.pressure_head_m.dtype, np.dtype("float64"))
        self.assertEqual(state.active.dtype, np.dtype(bool))
        self.assertEqual(tree[STATIC_GROUP].facies_id.dtype, np.dtype("int16"))
        self.assertNotIn("spatial_ref", state)
        self.assertNotIn("proj:wkt2", state.attrs)
        self.assertEqual(state.time.attrs["timezone"], "UTC")
        self.assertEqual(state.water_content.attrs["cell_methods"], "time: point")
        self.assertEqual(tree[EVENT_GROUP].event_time.values[0], np.datetime64("2026-01-01T00:00:37.123456789", "ns"))
        self.assertNotIn("time", tree[EVENT_GROUP].dims)
        np.testing.assert_array_equal(tree[INTERVAL_GROUP].time_bounds.values,
            np.array([["2026-01-01T00:00:00", "2026-01-01T01:00:00"]], dtype="datetime64[ns]"))

    def test_real_projected_crs_has_spatial_ref(self):
        tree = build_tree(**self.arguments(projected=True))
        self.assertEqual(tree[STATE_GROUP].to_dataset().rio.crs.to_epsg(), 32610)
        self.assertEqual(tree[STATE_GROUP].z.attrs["vertical_datum"], "NAVD88")

    def test_wrong_known_unit_rejected(self):
        args = self.arguments()
        args["state_fields"]["pressure_head_m"] = Field(np.zeros((2, 2, 2, 2)), "Pa")
        with self.assertRaisesRegex(ValueError, "units"):
            build_tree(**args)

    def test_missing_active_state_or_material_rejected(self):
        for group, name in (("state_fields", "pressure_head_m"), ("material_fields", "theta_s")):
            args = self.arguments()
            args[group][name].values[0, 0, 0, 0] = np.nan
            with self.subTest(group=group), self.assertRaisesRegex(ValueError, "finite.*active"):
                build_tree(**args)

    def test_physical_air_values_rejected(self):
        args = self.arguments()
        args["state_fields"]["pressure_head_m"].values[:, 1, 1, 1] = 0.
        with self.assertRaisesRegex(ValueError, "NaN outside"):
            build_tree(**args)

    def test_no_fabricated_standard_name(self):
        args = self.arguments()
        args["state_fields"]["pressure_head_m"] = Field(np.zeros((2, 2, 2, 2)), "m", {"standard_name": "fake_suction"})
        with self.assertRaisesRegex(ValueError, "standard_name"):
            build_tree(**args)

    def test_timezone_required_and_converted(self):
        args = self.arguments()
        args["state_times_utc"] = ["2025-12-31T16:00:00-08:00", "2025-12-31T17:00:00-08:00"]
        tree = build_tree(**args)
        self.assertEqual(tree[STATE_GROUP].time.values[0], np.datetime64("2026-01-01T00:00:00", "ns"))
        args["state_times_utc"] = ["2026-01-01T00:00:00", "2026-01-01T01:00:00"]
        with self.assertRaisesRegex(ValueError, "timezone"):
            build_tree(**args)
        args["state_times_utc"] = [np.datetime64("2026-01-01"), np.datetime64("2026-01-02")]
        with self.assertRaisesRegex(ValueError, "timezone"):
            build_tree(**args)

    def test_duplicate_snapshot_times_rejected(self):
        args = self.arguments()
        args["state_times_utc"] = ["2026-01-01T00:00:00Z"] * 2
        with self.assertRaisesRegex(ValueError, "unique"):
            build_tree(**args)

    def test_local_status_cannot_claim_real_crs(self):
        args = self.arguments(projected=True)
        args["reference"]["synthetic_local"] = True
        with self.assertRaisesRegex(ValueError, "disagree"):
            build_tree(**args)

    def test_synthetic_flag_requires_actual_boolean(self):
        for value in (0, 1, "true", None):
            args = self.arguments()
            args["reference"]["synthetic_local"] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "actual boolean"):
                canonical_reference(args["reference"])
        ref = self.arguments()["reference"]
        normalized = canonical_reference(ref)
        self.assertIs(normalized["synthetic_local"], True)
        self.assertNotIn("synthetic_local", ref)

    def test_metric_solver_rejects_geographic_or_feet_crs(self):
        for crs in ("EPSG:4326", "EPSG:2276"):
            args = self.arguments(projected=True)
            args["reference"]["horizontal_crs"] = crs
            with self.subTest(crs=crs), self.assertRaisesRegex(ValueError, "projected metre"):
                build_tree(**args)

    def test_transposed_state_shape_rejected(self):
        args = self.arguments()
        args["state_fields"]["pressure_head_m"] = Field(np.zeros((2, 2, 2)), "m")
        with self.assertRaisesRegex(ValueError, "shape"):
            build_tree(**args)

    def test_no_object_null_or_silent_zero_fill(self):
        args = self.arguments()
        args["static_fields"]["label"] = Field(np.full((2, 2, 2), None), "1")
        with self.assertRaisesRegex(ValueError, "only text"):
            build_tree(**args)

    def test_interval_requires_mean_or_sum_and_bounds(self):
        args = self.arguments()
        args["interval_fields"]["source_rate_m3_s"] = Field(np.array([1.]), "m3/s")
        with self.assertRaisesRegex(ValueError, "cell_methods"):
            build_tree(**args)
        args = self.arguments()
        args["interval_bounds_utc"] = None
        with self.assertRaisesRegex(ValueError, "explicit bounds"):
            build_tree(**args)

    def test_overlapping_intervals_rejected(self):
        args = self.arguments()
        args["interval_bounds_utc"] *= 2
        with self.assertRaisesRegex(ValueError, "overlap"):
            build_tree(**args)

    def test_source_license_and_source_link_required(self):
        args = self.arguments()
        args["provenance"]["links"] = [{"rel": "via", "href": "https://example.org"}]
        with self.assertRaisesRegex(ValueError, "license"):
            build_tree(**args)

    def test_state_material_epoch_mismatch_rejected(self):
        tree = build_tree(**self.arguments())
        ds = tree[MATERIAL_GROUP].to_dataset().assign_coords(material_epoch=("time", ["wrong", "wrong"]))
        tree[MATERIAL_GROUP] = ds
        with self.assertRaisesRegex(ValueError, "material_epoch"):
            validate_tree(tree)

    def test_material_mutation_cannot_keep_old_fingerprint(self):
        tree = build_tree(**self.arguments())
        tree[MATERIAL_GROUP]["theta_s"].values[1, 0, 0, 0] = .49
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            validate_tree(tree)

    def test_grid_mutation_cannot_keep_old_fingerprint(self):
        tree = build_tree(**self.arguments())
        state = tree[STATE_GROUP].to_dataset().assign_coords(x=[10., 11.])
        tree[STATE_GROUP] = state
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            validate_tree(tree)

    def test_event_unique_ids_and_order(self):
        args = self.arguments()
        args["event_ids"] = ["id", "id"]
        args["event_times_utc"] *= 2
        with self.assertRaisesRegex(ValueError, "unique"):
            build_tree(**args)

    def test_checkpoint_complete_memory_and_valid_time(self):
        args = self.arguments()
        tree = build_tree(**args)
        self.assertEqual(checkpoint_from_tree(tree), args["checkpoint_payload"])
        args["checkpoint_payload"]["valid_time_utc"] = "2026-01-01T00:00:00Z"
        with self.assertRaisesRegex(ValueError, "last model state"):
            build_tree(**args)

    def test_nonfinite_checkpoint_not_serialized(self):
        args = self.arguments()
        args["checkpoint_payload"]["bad"] = float("nan")
        with self.assertRaisesRegex(ValueError, "finite"):
            build_tree(**args)


def _gaia_available():
    try:
        from gaia_cli import io
        return io is not None
    except ImportError:
        return False


@unittest.skipIf(xr is None, "Install Hydrostate datahub dependencies")
class GaiaRoundtripTests(unittest.TestCase):
    arguments = DatahubTests.arguments

    def test_actual_gaia_typed_roundtrip_and_stac(self):
        if not _gaia_available():
            if os.environ.get("HYDROSTATE_REQUIRE_GAIA") == "1":
                self.fail("HYDROSTATE_REQUIRE_GAIA=1 requires actual Gaia io dependencies")
            self.skipTest("Install pinned gaia-cli to exercise its actual writer/reader")
        for projected in (False, True):
            with self.subTest(projected=projected), tempfile.TemporaryDirectory() as folder:
                tree = build_tree(**self.arguments(projected=projected))
                path = Path(folder) / "coupled.zarr"
                save_tree(tree, str(path))
                reopened = open_tree(str(path))
                state = reopened[STATE_GROUP].to_dataset()
                self.assertEqual(state.water_content.dtype, np.dtype("float32"))
                self.assertEqual(state.pressure_head_m.dtype, np.dtype("float64"))
                self.assertEqual(state.active.dtype, np.dtype(bool))
                self.assertFalse(bool(state.active.values[1, 1, 1]))
                self.assertTrue(np.isnan(state.water_content.values[:, 1, 1, 1]).all())
                self.assertEqual(reopened[STATIC_GROUP].facies_id.dtype, np.dtype("int16"))
                self.assertEqual(int(reopened[STATIC_GROUP].facies_id.values[1, 1, 1]), -1)
                np.testing.assert_array_equal(reopened[STATIC_GROUP].facies_label, tree[STATIC_GROUP].facies_label)
                np.testing.assert_array_equal(reopened[EVENT_GROUP].event_time, tree[EVENT_GROUP].event_time)
                np.testing.assert_array_equal(reopened[EVENT_GROUP].pressure_head_before_m, tree[EVENT_GROUP].pressure_head_before_m)
                self.assertEqual(checkpoint_from_tree(reopened), checkpoint_from_tree(tree))
                self.assertEqual(reopened.attrs["gaia:cli-version"], GAIA_COMMIT)
                metadata = json.loads((path / "zarr.json").read_text())
                self.assertEqual(metadata["zarr_format"], 3)
                self.assertIn("consolidated_metadata", metadata)
                product = reopened[STATE_GROUP].attrs
                self.assertEqual("proj:wkt2" in product, projected)
                item = json.loads((Path(folder) / "coupled.json").read_text())
                self.assertEqual(item["properties"]["license"], "MIT")
                self.assertEqual(item["geometry"] is not None, projected)
                self.assertEqual(item["properties"]["start_datetime"], "2026-01-01T00:00:00.000000000Z")
                reopened.close()


if __name__ == "__main__":
    unittest.main()
