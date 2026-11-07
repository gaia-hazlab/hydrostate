"""Analytical operator checks and CSV QA, independent of the flow integrator."""

from __future__ import annotations

from copy import deepcopy
import csv
from dataclasses import replace
from datetime import datetime, timedelta, timezone
import json
import math
from pathlib import Path
import tempfile
import unittest

from hydrostate.observations import (
    CANONICAL_UNITS, OBSERVATION_COLUMNS, SENSOR_COLUMNS, ModelSampler,
    Observation, ObservationError, Sensor, compare_observations, iso_utc,
    load_observations, load_sensors, main, metrics, unit_factor, utc_time,
)


START = datetime(2026, 1, 1, tzinfo=timezone.utc)
CRS = "synthetic local Cartesian metres"
DATUM = "arbitrary local elevation datum"


def analytical_model() -> dict:
    axis = [0.5, 1.5, 2.5]
    snapshots = []
    for time in (0, 1, 2):
        def field(function):
            return [[[function(x, y, z, time) for x in axis] for y in axis] for z in axis]
        snapshots.append({
            "time_h": time,
            "theta": field(lambda x, y, z, t: 0.1 + 0.01 * x + 0.02 * y + 0.03 * z + 0.004 * t),
            "pressure_head_m": field(lambda x, y, z, t: 10 - z + 0.1 * x + 0.2 * y + 0.5 * t),
            "total_head_m": field(lambda x, y, z, t: 10 + 0.1 * x + 0.2 * y + 0.5 * t),
        })
    return {
        "meta": {"start_utc": iso_utc(START), "coordinate_reference": {"horizontal_crs": CRS, "vertical_datum": DATUM}, "aquifer_id": "one_aquifer", "precipitation_spatial_support": "uniform", "outlet": {"id": "east", "x_m": 3.0, "y_m": 1.5}},
        "grid": {"x_m": axis, "y_m": axis, "z_m": axis, "dx_m": 1.0, "dy_m": 1.0, "dz_m": 1.0, "active": [[[True] * 3 for _ in axis] for _ in axis], "surface_elevation_m": [[3.0] * 3 for _ in axis]},
        "snapshots": snapshots,
        "records": [{"time_h": 0, "outlet_m3_s": 999, "precip_mm_h": 999}, {"time_h": 1, "outlet_m3_s": 2, "precip_mm_h": 10}, {"time_h": 2, "outlet_m3_s": 6, "precip_mm_h": 30}],
    }


def sensor(kind="soil_moisture", **kwargs) -> Sensor:
    values = dict(sensor_id="S", kind=kind, x_m=1.0, y_m=1.25, z_m=1.75,
                  units=CANONICAL_UNITS[kind], uncertainty=0.1,
                  horizontal_crs=CRS, vertical_datum=DATUM, role="validation",
                  provenance="synthetic analytical test")
    values.update(kwargs)
    return Sensor(**values)


def observation(s: Sensor, time=1.0, value=0.2, start=None, **kwargs) -> Observation:
    values = dict(sensor_id=s.sensor_id, valid_time=START + timedelta(hours=time),
                  arrival_time=START + timedelta(hours=time, minutes=10),
                  value=value, uncertainty=0.1, units=CANONICAL_UNITS[s.kind],
                  quality_flag="good", interval_start=START + timedelta(hours=start) if start is not None else None)
    values.update(kwargs)
    return Observation(**values)


class OperatorTests(unittest.TestCase):
    def setUp(self):
        self.model = analytical_model()
        self.sampler = ModelSampler(self.model)

    def test_linear_field_point_and_time_interpolation(self):
        s = sensor()
        actual = self.sampler.predict(s, observation(s, time=0.5))
        expected = 0.1 + 0.01 * 1.0 + 0.02 * 1.25 + 0.03 * 1.75 + 0.004 * 0.5
        self.assertAlmostEqual(actual, expected, places=13)

    def test_support_average_and_depth_below_surface(self):
        s = sensor(x_m=1.5, y_m=1.5, z_m=None, depth_m=1.5,
                   support_x_m=2, support_y_m=2, support_z_m=2)
        self.assertAlmostEqual(self.sampler.predict(s, observation(s)), 0.1 + 0.06 * 1.5 + 0.004)

    def test_pore_pressure_is_pressure_head_times_rho_g(self):
        s = sensor("pore_pressure")
        psi = 10 - 1.75 + 0.1 * 1 + 0.2 * 1.25 + 0.5
        self.assertAlmostEqual(self.sampler.predict(s, observation(s)), psi * 9806.65)

    def test_hydraulic_head_is_not_pressure_head(self):
        s = sensor("hydraulic_head")
        head = self.sampler.predict(s, observation(s))
        psi_sensor = replace(s, kind="pressure_head")
        psi = self.sampler.predict(psi_sensor, observation(psi_sensor))
        self.assertAlmostEqual(head - psi, 1.75)

    def test_screen_length_weights_and_level_is_not_water_table(self):
        for frame in self.model["snapshots"]:
            frame["total_head_m"] = [[[10 + z for _ in range(3)] for _ in range(3)] for z in [0.5, 1.5, 2.5]]
        s = sensor("hydraulic_head", x_m=1.5, y_m=1.5, z_m=None,
                   screen_bottom_m=0.25, screen_top_m=2.25, aquifer_id="one_aquifer",
                   screen_assumption="single_aquifer_uniform_conductance")
        head = ModelSampler(self.model).predict(s, observation(s))
        self.assertAlmostEqual(head, 11.25)
        level = replace(s, kind="water_level_depth")
        self.assertAlmostEqual(ModelSampler(self.model).predict(level, observation(level)), 3 - head)
        self.assertLess(3 - head, 0)  # Artesian piezometric level is not clipped to zero.

    def test_unsaturated_or_wrong_aquifer_screen_rejected(self):
        s = sensor("hydraulic_head", x_m=0.5, y_m=0.5, z_m=None,
                   screen_bottom_m=0, screen_top_m=1, aquifer_id="wrong",
                   screen_assumption="single_aquifer_uniform_conductance")
        with self.assertRaisesRegex(ObservationError, "aquifer"):
            self.sampler.predict(s, observation(s))
        s = replace(s, aquifer_id="one_aquifer")
        self.model["snapshots"][1]["pressure_head_m"][0][0][0] = -0.1
        with self.assertRaisesRegex(ObservationError, "unsaturated"):
            ModelSampler(self.model).predict(s, observation(s))

    def test_rate_operator_integrates_preceding_hour_means(self):
        s = sensor("streamflow", z_m=None, outlet_id="east")
        self.assertAlmostEqual(self.sampler.predict(s, observation(s, time=1.5, start=0.5)), 4.0)
        self.assertAlmostEqual(self.sampler.predict(s, observation(s, time=1, start=0)), 2.0)
        self.assertAlmostEqual(self.sampler.predict(s, observation(s, time=2, start=0)), 4.0)

    def test_precipitation_point_requires_uniform_forcing(self):
        s = sensor("precipitation", z_m=None)
        self.assertAlmostEqual(self.sampler.predict(s, observation(s, time=1.75, start=0.25)), 20.0)
        del self.model["meta"]["precipitation_spatial_support"]
        with self.assertRaisesRegex(ObservationError, "point-gauge"):
            ModelSampler(self.model).predict(s, observation(s, time=1, start=0))

    def test_mismatched_reference_outside_grid_or_time_rejected(self):
        cases = [(sensor(vertical_datum="NAVD88"), 1, "datum"),
                 (sensor(horizontal_crs="EPSG:32610"), 1, "CRS"),
                 (sensor(x_m=3.1), 1, "outside"),
                 (sensor(), 3, "outside")]
        for s, time, pattern in cases:
            with self.subTest(pattern=pattern), self.assertRaisesRegex(ObservationError, pattern):
                self.sampler.predict(s, observation(s, time=time))

    def test_inactive_support_not_renormalized(self):
        self.model["grid"]["active"][1][1][1] = False
        with self.assertRaisesRegex(ObservationError, "inactive"):
            ModelSampler(self.model).predict(sensor(), observation(sensor()))

    def test_role_and_forcing_metrics_remain_separate(self):
        a = sensor(sensor_id="A", role="train")
        b = replace(a, sensor_id="B", role="validation")
        p = sensor("precipitation", sensor_id="P", z_m=None)
        obs = [observation(a), observation(b), observation(p, start=0, value=10)]
        report = compare_observations(self.model, {s.sensor_id: s for s in (a, b, p)}, obs)
        self.assertEqual({r["role"] for r in report["state_and_discharge_validation"]}, {"train", "validation"})
        self.assertEqual([r["kind"] for r in report["forcing_validation"]], ["precipitation"])
        self.assertIn("no assimilation", report["method"])
        self.assertIn("not independent field validation", report["scope"])

    def test_missing_qc_and_arrival_cutoff_are_exclusions_not_zeros(self):
        s = sensor()
        obs = [observation(s, value=None), observation(s, time=2, quality_flag="bad")]
        report = compare_observations(self.model, {"S": s}, obs)
        self.assertEqual(report["matched_count"], 0)
        self.assertEqual({v["reason"] for v in report["excluded"]}, {"missing_value", "quality_flag"})
        report = compare_observations(self.model, {"S": s}, [observation(s)], as_of=START + timedelta(hours=1, minutes=5))
        self.assertEqual(report["excluded"][0]["reason"], "not_available_as_of")

    def test_weighted_metrics_detect_known_perturbations(self):
        rows = [{"predicted": 5, "observed": 4, "uncertainty": 1}, {"predicted": 2, "observed": 4, "uncertainty": 2}]
        result = metrics(rows)
        self.assertAlmostEqual(result["bias_model_minus_observed"], -0.5)
        self.assertAlmostEqual(result["rmse"], math.sqrt(2.5))
        self.assertAlmostEqual(result["uncertainty_weighted_rmse"], math.sqrt(1.6))
        self.assertAlmostEqual(result["mean_squared_observation_normalized_residual"], 1.0)


class CsvTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base = Path(self.directory.name)

    def write(self, name, columns, rows):
        path = self.base / name
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader(); writer.writerows(rows)
        return path

    def metadata(self, **changes):
        row = {"sensor_id": "S", "kind": "soil_moisture", "x_m": 1.0, "y_m": 1.25,
               "z_m": 1.75, "units": "%", "uncertainty": 2, "role": "validation",
               "horizontal_crs": CRS, "vertical_datum": DATUM, "provenance": "synthetic"}
        row.update(changes)
        return row

    def measurement(self, **changes):
        row = {"sensor_id": "S", "valid_time_utc": "2026-01-01T01:00:00Z", "arrival_time_utc": "2026-01-01T01:10:00Z", "value": 30, "units": "%", "quality_flag": "good"}
        row.update(changes)
        return row

    def test_percent_and_uncertainty_normalization(self):
        sensors = load_sensors(self.write("s.csv", SENSOR_COLUMNS, [self.metadata()]))
        obs = load_observations(self.write("o.csv", OBSERVATION_COLUMNS, [self.measurement()]), sensors)
        self.assertEqual(obs[0].value, 0.3)
        self.assertEqual(obs[0].uncertainty, 0.02)
        self.assertEqual(obs[0].units, "m3/m3")

    def test_duplicate_metadata_and_normalized_timestamps_rejected(self):
        with self.assertRaisesRegex(ObservationError, "duplicate"):
            load_sensors(self.write("s.csv", SENSOR_COLUMNS, [self.metadata(), self.metadata()]))
        sensors = load_sensors(self.write("s.csv", SENSOR_COLUMNS, [self.metadata()]))
        rows = [self.measurement(), self.measurement(valid_time_utc="2025-12-31T17:00:00-08:00")]
        with self.assertRaisesRegex(ObservationError, "Duplicate"):
            load_observations(self.write("o.csv", OBSERVATION_COLUMNS, rows), sensors)

    def test_naive_time_bad_quality_nonfinite_and_bounds_rejected(self):
        sensors = load_sensors(self.write("s.csv", SENSOR_COLUMNS, [self.metadata()]))
        for changes in [{"valid_time_utc": "2026-01-01T01:00:00"}, {"quality_flag": "probably"}, {"value": "inf"}, {"value": 120}, {"arrival_time_utc": "2026-01-01T00:00:00Z"}]:
            with self.subTest(changes=changes), self.assertRaises(ObservationError):
                load_observations(self.write("o.csv", OBSERVATION_COLUMNS, [self.measurement(**changes)]), sensors)

    def test_blank_value_remains_missing(self):
        sensors = load_sensors(self.write("s.csv", SENSOR_COLUMNS, [self.metadata()]))
        obs = load_observations(self.write("o.csv", OBSERVATION_COLUMNS, [self.measurement(value="", quality_flag="missing")]), sensors)
        self.assertIsNone(obs[0].value)

    def test_duplicate_header_and_missing_required_columns_rejected(self):
        path = self.base / "bad.csv"
        path.write_text("sensor_id,sensor_id\nA,B\n")
        with self.assertRaisesRegex(ObservationError, "duplicate CSV header"):
            load_sensors(path)
        path.write_text("sensor_id,kind\nA,soil_moisture\n")
        with self.assertRaisesRegex(ObservationError, "missing columns"):
            load_sensors(path)

    def test_precip_accumulation_interval_and_conversion(self):
        sensors = load_sensors(self.write("s.csv", SENSOR_COLUMNS, [self.metadata(kind="precipitation", z_m="", units="mm", uncertainty=1)]))
        row = self.measurement(value=6, units="mm", interval_start_utc="2026-01-01T00:30:00Z")
        obs = load_observations(self.write("o.csv", OBSERVATION_COLUMNS, [row]), sensors)[0]
        self.assertEqual(obs.value, 12)
        self.assertEqual(obs.uncertainty, 2)
        self.assertAlmostEqual(unit_factor("streamflow", "cfs"), 0.028316846592)
        self.assertEqual(unit_factor("hydraulic_head", "ft"), 0.3048)

    def test_cli_demo_and_validation_round_trip(self):
        path = self.base / "model.json"
        path.write_text(json.dumps(analytical_model()))
        output = self.base / "report.json"
        self.assertEqual(main(["--model", str(path), "--make-demo", str(self.base / "demo"), "--output", str(output)]), 0)
        report = json.loads(output.read_text())
        self.assertGreater(report["matched_count"], 0)
        self.assertEqual(report["excluded_count"], 1)
        self.assertGreater(next(r for r in report["by_sensor"] if r["sensor_id"] == "SM_2")["rmse"], 0)

    def test_utc_offset_normalization(self):
        self.assertEqual(utc_time("2025-12-31T17:00:00-08:00"), START + timedelta(hours=1))


if __name__ == "__main__":
    unittest.main()
