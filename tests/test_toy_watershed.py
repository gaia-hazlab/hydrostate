"""Numerical reference checks, not field calibration or hazard validation."""

import json
import importlib.util
import unittest

if importlib.util.find_spec("numpy") is not None and importlib.util.find_spec("scipy") is not None:
    import numpy as np
    from puget_modflow.toy_watershed import RichardsHillslope, ToyConfig, retention, run_toy
else:
    np = None


@unittest.skipIf(np is None, "toy reference requires optional numpy and scipy")
class ToyConstitutiveTests(unittest.TestCase):
    def test_retention_bounds_monotonicity_and_capacity(self):
        psi = np.r_[-np.logspace(5, -5, 150), 0.0, 0.3, 3.0]
        theta, capacity, kr = retention(psi)
        self.assertTrue(np.all(theta >= 0.06))
        self.assertTrue(np.all(theta <= 0.43))
        self.assertTrue(np.all(np.diff(theta) >= 0))
        self.assertTrue(np.all((kr >= 0) & (kr <= 1)))
        self.assertTrue(np.all(capacity >= 0))
        self.assertTrue(np.all(theta[psi >= 0] == 0.43))
        sample = np.array([-20.0, -2.0, -0.1])
        epsilon = 1.0e-6
        derivative = (retention(sample + epsilon)[0] - retention(sample - epsilon)[0]) / (2 * epsilon)
        np.testing.assert_allclose(derivative, retention(sample)[1], rtol=1.0e-6)

    def test_rotated_hydrostatic_no_flow(self):
        model = RichardsHillslope(ToyConfig(nx=6, nz=12, hours=1))
        psi = model.initial_pressure()
        surface = np.zeros(model.config.nx)
        result, surf, stage, info = model.advance(
            psi, surface, model.initial_head, 3600, closed=True
        )
        np.testing.assert_allclose(result, psi, atol=1.0e-12)
        np.testing.assert_allclose(result + model.z, model.initial_head, atol=1.0e-12)
        np.testing.assert_array_equal(surf, surface)
        self.assertEqual(stage, model.initial_head)
        self.assertLess(info["local_max_residual_m3_s"], 1.0e-12)

    def test_closed_column_capillary_redistribution_conserves_water(self):
        config = ToyConfig(nx=4, nz=12, length_m=20, depth_m=3, width_m=5, slope_rise_run=0, hours=1)
        model = RichardsHillslope(config, lens=False)
        # Uniform suction is NOT hydrostatic: z+psi varies with depth. With all
        # faces sealed, gravity drains the top and wets the bottom, conserving
        # the entire column's water as suction adjusts.
        old = np.full(model.z.shape, -1.0)
        before = float(model.storage_density(old).sum() * model.volume)
        psi, _, _, info = model.advance(old, np.zeros(config.nx), model.initial_head, 600, closed=True)
        after = float(model.storage_density(psi).sum() * model.volume)
        self.assertLess(psi[0, 0], old[0, 0])
        self.assertGreater(psi[-1, 0], old[-1, 0])
        self.assertLess(abs(after - before), 1.0e-6)
        self.assertLess(info["final_pressure_increment_m"], config.pressure_tolerance_m)

    def test_wet_channel_loses_water_to_lower_head_soil(self):
        config = ToyConfig(nx=4, nz=12, length_m=20, depth_m=3, width_m=5, slope_rise_run=0, hours=1)
        model = RichardsHillslope(config, lens=False)
        old = model.initial_pressure()
        stage_old = model.initial_head + 0.5
        soil_old = float(model.storage_density(old).sum() * model.volume)
        channel_old = config.channel_area_m2 * (stage_old - model.channel_bed)
        dt = 60.0
        psi, surface, stage, info = model.advance(old, np.zeros(config.nx), stage_old, dt)
        soil_new = float(model.storage_density(psi).sum() * model.volume)
        channel_new = config.channel_area_m2 * (stage - model.channel_bed)
        outlet = channel_new / config.channel_tau_s * dt
        self.assertLess(info["qtoe_m3_s"], 0)
        self.assertGreater(soil_new, soil_old)
        residual = soil_new - soil_old + channel_new - channel_old + outlet + surface.sum() * model.surface_area
        self.assertLess(abs(residual), 1.0e-6)

    def test_channel_does_not_exchange_through_below_bed_wall(self):
        model = RichardsHillslope(ToyConfig(nx=4, nz=24, hours=1))
        psi = model.initial_pressure()
        # Pressure in the deep closed toe wall may affect internal gradients,
        # but cannot leak directly into an above-bed bank reservoir.
        wholly_below_bed = model.toe_z + 0.5 * model.dz * model.cos_beta < model.channel_bed
        self.assertTrue(np.any(wholly_below_bed))
        psi[wholly_below_bed, -1] += 1.0
        state = np.r_[psi.ravel(), np.zeros(model.config.nx), model.initial_head]
        _, _, info = model._system(state, state, 60.0, 0.0, 0.0)
        self.assertAlmostEqual(info["qtoe_m3_s"], 0.0, places=14)


@unittest.skipIf(np is None, "toy reference requires optional numpy and scipy")
class ToyStormReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Includes snowfall, warm rain-on-snow, a later intense pulse, top
        # ponding/runoff, and both directions of wet bank exchange.
        cls.coarse = run_toy(hours=72, snapshot_interval_h=6)
        cls.fine = run_toy(hours=72, max_step_s=300, snapshot_interval_h=6)

    def test_all_stores_close_and_local_residual_is_small(self):
        result = self.coarse
        self.assertTrue(result["budget"]["accepted"])
        self.assertLess(abs(result["budget"]["signed_residual_mm"]), 1.0e-5)
        self.assertLess(result["budget"]["max_abs_hourly_cumulative_residual_m3"], 0.001)
        self.assertLess(result["diagnostics"]["max_local_mass_residual_m3_s"], 2.0e-9)
        self.assertLess(result["diagnostics"]["max_accepted_pressure_increment_m"], 1.0e-5)
        self.assertGreater(max(row["runoff_m3_s"] for row in result["records"]), 0)
        self.assertGreater(max(row["melt_mm_h"] for row in result["records"]), 0)
        self.assertGreater(result["records"][-1]["theta_root"], result["records"][0]["theta_root"])
        # Reject nonfinite output instead of silently producing nonstandard JSON.
        json.dumps(result, allow_nan=False)

    def test_hourly_states_are_insensitive_to_halved_max_step(self):
        a, b = self.coarse["records"], self.fine["records"]
        self.assertEqual(len(a), len(b))
        max_difference = lambda key: max(abs(x[key] - y[key]) for x, y in zip(a, b))
        self.assertLess(max_difference("theta_root"), 0.001)
        self.assertLess(max_difference("wt_depth_m"), 0.015)
        self.assertLess(max_difference("pressure_sensor_kpa"), 0.15)
        self.assertLess(max_difference("outlet_m3_s"), 0.0002)

    def test_water_table_diagnostic_excludes_an_isolated_perch(self):
        model = RichardsHillslope(ToyConfig(nx=4, nz=12, hours=1))
        psi = np.full(model.z.shape, -1.0)
        psi[3:5] = 0.2
        self.assertEqual(model.water_table_depth(psi), [None] * 4)
        psi[-3:] = 0.3
        table = model.water_table_depth(psi)
        self.assertTrue(all(d is not None and d > 4 for d in table))


if __name__ == "__main__":
    unittest.main()
