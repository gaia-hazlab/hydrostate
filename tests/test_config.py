from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from puget_modflow.config import PipelineConfig, write_default_config


class ConfigTests(unittest.TestCase):
    def test_default_config_loads(self) -> None:
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "puget.toml"
            write_default_config(path)
            config = PipelineConfig.from_toml(path)

        self.assertEqual(config.domain.name, "southeast_puget_sound")
        self.assertEqual(config.model.nlay, 13)
        self.assertGreaterEqual(len(config.reference_datasets), 3)
        self.assertEqual(config.precipitation.element, "ppt")


if __name__ == "__main__":
    unittest.main()
