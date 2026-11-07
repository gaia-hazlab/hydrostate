import unittest

from puget_modflow.datasources import (
    month_tokens,
    nwis_groundwater_levels_url,
    prism_data_url,
    prism_release_date_url,
    sciencebase_item_json_url,
    sciencebase_attached_files,
)


class DataSourceTests(unittest.TestCase):
    def test_prism_monthly_url(self) -> None:
        url = prism_data_url("ppt", "200501", region="us", resolution="4km")
        self.assertEqual(url, "https://services.nacse.org/prism/data/get/us/4km/ppt/200501")

    def test_prism_url_with_format(self) -> None:
        url = prism_data_url("ppt", "2005-01-15", fmt="nc")
        self.assertEqual(url, "https://services.nacse.org/prism/data/get/us/4km/ppt/20050115?format=nc")

    def test_release_url(self) -> None:
        url = prism_release_date_url("ppt", "200501")
        self.assertIn("/releaseDate/us/4km/ppt/200501", url)
        self.assertTrue(url.endswith("json=true"))

    def test_month_tokens(self) -> None:
        self.assertEqual(month_tokens("2005-11", "2006-02"), ["200511", "200512", "200601", "200602"])

    def test_sciencebase_url(self) -> None:
        self.assertEqual(
            sciencebase_item_json_url("62fa7dfed34e25bd092879e6"),
            "https://www.sciencebase.gov/catalog/item/62fa7dfed34e25bd092879e6?format=json",
        )

    def test_nwis_url(self) -> None:
        url = nwis_groundwater_levels_url(
            bbox=(-123.3, 46.75, -121.7, 47.55),
            start="2005-01-01",
            end="2015-12-31",
        )
        self.assertIn("siteType=GW", url)
        self.assertIn("parameterCd=72019", url)
        self.assertIn("startDT=2005-01-01", url)
        self.assertIn("endDT=2015-12-31", url)

    def test_sciencebase_attached_files(self) -> None:
        item = {
            "files": [
                {"name": "model.zip", "downloadUri": "https://example.test/model.zip"},
                {"name": "", "downloadUri": "https://example.test/skip.zip"},
            ]
        }
        self.assertEqual(
            sciencebase_attached_files(item),
            [{"name": "model.zip", "url": "https://example.test/model.zip"}],
        )


if __name__ == "__main__":
    unittest.main()
