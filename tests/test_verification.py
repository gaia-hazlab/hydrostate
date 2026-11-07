"""Analytical verification gates for hydrostate's actual 3-D flow kernel."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

try:
    from hydrostate.verification import run_verification
except ModuleNotFoundError as error:
    if error.name not in {"numpy", "scipy"}:
        raise
    run_verification = None


@unittest.skipIf(run_verification is None, "3-D flow verification requires numpy and scipy")
class MathematicalVerificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = run_verification()

    def test_constitutive_storage_and_conductivity_derivatives(self):
        self.assertTrue(self.report["checks"]["constitutive_derivatives"]["passed"])

    def test_affine_darcy_fluxes_in_three_dimensions(self):
        check = self.report["checks"]["affine_darcy_head"]
        self.assertTrue(check["passed"], check)

    def test_source_and_prescribed_flux_units_and_sign(self):
        check = self.report["checks"]["prescribed_source_and_flux"]
        self.assertTrue(check["passed"], check)
        self.assertLess(check["storage_change_m3"], 0)

    def test_hydrostatic_equilibrium_on_heterogeneous_active_mask(self):
        check = self.report["checks"]["heterogeneous_masked_hydrostatic"]
        self.assertGreater(check["inactive_cells"], 0)
        self.assertTrue(check["passed"], check)

    def test_nonlinear_closed_domain_redistributes_and_conserves(self):
        check = self.report["checks"]["nonlinear_closed_redistribution"]
        self.assertTrue(check["passed"], check)

    def test_actual_linear_assembly_matches_residual_derivative(self):
        check = self.report["checks"]["linear_residual_jacobian"]
        self.assertTrue(check["passed"], check)

    def test_three_dimensional_mms_has_second_order_space(self):
        check = self.report["checks"]["linear_3d_spatial_mms"]
        self.assertTrue(check["passed"], check["rates"])

    def test_exact_time_source_gives_first_order_backward_euler(self):
        check = self.report["checks"]["backward_euler_temporal_mms"]
        self.assertTrue(check["passed"], check["rates"])

    def test_smooth_nonlinear_mms_remains_an_explicit_diagnostic(self):
        check = self.report["checks"]["nonlinear_3d_spatial_mms"]
        self.assertIsNone(check["passed"])
        self.assertFalse(check["acceptance_required"])
        self.assertFalse(check["asymptotic_linf_verified"])
        self.assertTrue(all(b["l2_m"] < a["l2_m"]
                            for a, b in zip(check["records"], check["records"][1:])))

    def test_complete_report_is_finite_identifies_source_and_passes(self):
        json.dumps(self.report, allow_nan=False)
        self.assertTrue(self.report["passed"])
        self.assertEqual(len(self.report["kernel_sha256"]), 64)
        self.assertIn("not field validation", self.report["limitations"][0])

    def test_failed_required_acceptance_returns_nonzero_cli_status(self):
        from contextlib import redirect_stdout
        import io
        from hydrostate.verification import main
        failure = {"passed": False, "checks": {"sentinel": {"passed": False}}, "wall_time_s": 0.0}
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "failed.json"
            with patch("hydrostate.verification.run_verification", return_value=failure):
                with redirect_stdout(io.StringIO()):
                    status = main(["--output", str(output)])
            self.assertEqual(status, 1)
            self.assertFalse(json.loads(output.read_text())["passed"])


if __name__ == "__main__":
    unittest.main()
