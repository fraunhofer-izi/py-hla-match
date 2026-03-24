import unittest
import os
from py_hla_match.external import query_dpb1_tce, DPB1TCEStatus


class TestDPB1TCEIntegration(unittest.TestCase):
    """
    Integration tests against the real EBI DPB1 TCE API.

    These tests are run weekly to ensure the external API is still functional.
    They use real API calls with known allele combinations.
    """

    def setUp(self):
        """Skip tests if running in CI without integration flag."""
        if not os.getenv('RUN_INTEGRATION_TESTS'):
            self.skipTest("Integration tests disabled")

    def test_api_v21_permissive_match(self):
        """Test API v2.1 with known permissive alleles."""
        result = query_dpb1_tce(
            patient_dpb1="01:01",
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="2.1",
            timeout=30
        )

        self.assertEqual(
            result.status,
            DPB1TCEStatus.SUCCESS,
            f"API v2.1 failed with status: {result.status}"
        )
        self.assertIsNotNone(result.prediction)
        self.assertTrue(
            "Permissive" in result.prediction or
            "ARD Matched" in result.prediction,
            f"Expected Permissive or ARD Matched, got: '{result.prediction}'"
        )

    def test_api_v21_non_permissive_gvh(self):
        """Test API v2.1 with known non-permissive alleles."""
        result = query_dpb1_tce(
            patient_dpb1="01:01",
            patient_dpb2="02:01",
            donor_dpb1="03:01",
            donor_dpb2="04:01",
            version="2.1",
            timeout=30
        )

        self.assertEqual(
            result.status,
            DPB1TCEStatus.SUCCESS,
            f"API v2.1 failed with status: {result.status}"
        )

        self.assertIsNotNone(result.prediction)

        prediction_normalized = result.prediction.lower().replace("-", " ")

        self.assertIn(
            "non permissive",
            prediction_normalized,
            f"Expected Non-Permissive prediction, got: '{result.prediction}'"
        )

    def test_api_v20_still_functional(self):
        """Test that older API v2.0 still works."""
        result = query_dpb1_tce(
            patient_dpb1="01:01",
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="2.0",
            timeout=30
        )

        self.assertEqual(
            result.status,
            DPB1TCEStatus.SUCCESS,
            f"API v2.0 failed with status: {result.status}"
        )

    def test_api_v30_still_functional(self):
        """Test that newer API v3.0 still works."""
        result = query_dpb1_tce(
            patient_dpb1="01:01",
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="3.0",
            timeout=30
        )

        self.assertEqual(
            result.status,
            DPB1TCEStatus.SUCCESS,
            f"API v3.0 failed with status: {result.status}"
        )
        self.assertIsNotNone(result.prediction)

    def test_invalid_alleles_handling(self):
        """Test that API properly rejects invalid alleles with HTTP 422."""
        result = query_dpb1_tce(
            patient_dpb1="99:99",  # Invalid allele
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="2.1",
            timeout=30
        )

        # Returns HTTP 422 which becomes REQUEST_ERROR in our wrapper
        self.assertEqual(
            result.status,
            DPB1TCEStatus.REQUEST_ERROR,
            f"Expected REQUEST_ERROR (422), got {result.status}"
        )


if __name__ == '__main__':
    # Allow running locally with:
    # RUN_INTEGRATION_TESTS=1 python -m pytest tests/test_integration.py
    os.environ['RUN_INTEGRATION_TESTS'] = '1'
    unittest.main()
