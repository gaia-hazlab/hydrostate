from datetime import date
import unittest

from puget_modflow.precip import precip_mm_to_inches, prism_date_from_name


class PrecipTests(unittest.TestCase):
    def test_mm_to_inches(self) -> None:
        self.assertAlmostEqual(precip_mm_to_inches(25.4), 1.0)

    def test_prism_monthly_date_from_name(self) -> None:
        self.assertEqual(
            prism_date_from_name("prism_ppt_us_4km_200501.zip"),
            date(2005, 1, 1),
        )

    def test_prism_daily_date_from_name(self) -> None:
        self.assertEqual(
            prism_date_from_name("PRISM_ppt_stable_4kmD2_20050115_bil.zip"),
            date(2005, 1, 15),
        )


if __name__ == "__main__":
    unittest.main()
