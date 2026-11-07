from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from puget_modflow.config import RechargeConfig
from puget_modflow.recharge import (
    MonthlyPrecipRecord,
    monthly_precip_to_recharge_rate_ft_per_day,
    read_monthly_precip_csv,
    recharge_series_from_monthly_precip,
    write_recharge_timeseries,
)


class RechargeTests(unittest.TestCase):
    def test_monthly_precip_to_daily_recharge(self) -> None:
        value = monthly_precip_to_recharge_rate_ft_per_day(
            precip_in=12.0,
            year=2005,
            month=1,
            fraction=0.25,
        )
        self.assertAlmostEqual(value, 0.25 / 31.0)

    def test_recharge_series_uses_seasonal_factor(self) -> None:
        config = RechargeConfig(
            method="test",
            default_fraction=0.2,
            min_fraction=0.0,
            max_fraction=1.0,
            units="feet_per_day",
            seasonal_factor={"jan": 2.0},
        )
        records = [MonthlyPrecipRecord(date=date(2005, 1, 1), precip_in=6.0)]
        recharge = recharge_series_from_monthly_precip(records, config)
        self.assertAlmostEqual(recharge[0].recharge_ft_per_day, (0.5 * 0.2 * 2.0) / 31.0)

    def test_csv_roundtrip(self) -> None:
        with TemporaryDirectory() as tmp:
            precip_path = Path(tmp) / "precip.csv"
            precip_path.write_text("date,precip_in\n2005-01,4.2\n", encoding="utf-8")
            precip = read_monthly_precip_csv(precip_path)
            config = RechargeConfig(
                method="test",
                default_fraction=0.2,
                min_fraction=0.0,
                max_fraction=1.0,
                units="feet_per_day",
                seasonal_factor={"jan": 1.0},
            )
            output = Path(tmp) / "recharge.csv"
            write_recharge_timeseries(recharge_series_from_monthly_precip(precip, config), output)
            text = output.read_text(encoding="utf-8")

        self.assertIn("recharge_ft_per_day", text)
        self.assertIn("2005-01-01", text)

    def test_csv_accepts_prism_mm(self) -> None:
        with TemporaryDirectory() as tmp:
            precip_path = Path(tmp) / "precip.csv"
            precip_path.write_text("date,precip_mm\n2005-01,25.4\n", encoding="utf-8")
            precip = read_monthly_precip_csv(precip_path)

        self.assertAlmostEqual(precip[0].precip_in, 1.0)


if __name__ == "__main__":
    unittest.main()
