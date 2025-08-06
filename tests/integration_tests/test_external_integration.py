# tests/test_integration.py
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
        """Skip tests if running in CI without integration flag"""
        if not os.getenv('RUN_INTEGRATION_TESTS'):
            self.skipTest("Integration tests disabled")

    def test_api_v21_permissive_match(self):
        """Test API v2.1 with known permissive alleles"""
        result = query_dpb1_tce(
            patient_dpb1="01:01",
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="2.1",
            timeout=30
        )

        # Should get result (not error)
        self.assertNotIn(result, [
            DPB1TCEStatus.API_ERROR,
            DPB1TCEStatus.TIMEOUT_ERROR,
            DPB1TCEStatus.REQUEST_ERROR,
            DPB1TCEStatus.VALUE_ERROR,
            DPB1TCEStatus.UNEXPECTED_ERROR
        ], f"API v2.1 returned error status: {result}")

        # For identical alleles, expect permissive
        self.assertEqual(
            result,
            DPB1TCEStatus.PERMISSIVE,
            "Identical alleles should be permissive"
        )

    def test_api_v21_non_permissive_gvh(self):
        """Test API v2.1 with known non-permissive GvH alleles"""
        result = query_dpb1_tce(
            patient_dpb1="02:01",
            patient_dpb2="02:01",
            donor_dpb1="04:02",
            donor_dpb2="04:02",
            version="2.1",
            timeout=30
        )

        # Should get result (not error)
        self.assertNotIn(result, [
            DPB1TCEStatus.API_ERROR,
            DPB1TCEStatus.TIMEOUT_ERROR,
            DPB1TCEStatus.REQUEST_ERROR,
            DPB1TCEStatus.VALUE_ERROR,
            DPB1TCEStatus.UNEXPECTED_ERROR
        ], f"API v2.1 returned error status: {result}")

    def test_api_v20_still_functional(self):
        """Test that older API v2.0 still works"""
        result = query_dpb1_tce(
            patient_dpb1="01:01",
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="2.0",
            timeout=30
        )

        # Should get result (not error)
        self.assertNotIn(result, [
            DPB1TCEStatus.API_ERROR,
            DPB1TCEStatus.TIMEOUT_ERROR,
            DPB1TCEStatus.REQUEST_ERROR,
            DPB1TCEStatus.VALUE_ERROR,
            DPB1TCEStatus.UNEXPECTED_ERROR
        ], f"API v2.0 returned error status: {result}")

    def test_api_v30_still_functional(self):
        """Test that newer API v3.0 still works"""
        result = query_dpb1_tce(
            patient_dpb1="01:01",
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="3.0",
            timeout=30
        )
        # Should get result (not error)
        self.assertNotIn(result, [
            DPB1TCEStatus.API_ERROR,
            DPB1TCEStatus.TIMEOUT_ERROR,
            DPB1TCEStatus.REQUEST_ERROR,
            DPB1TCEStatus.VALUE_ERROR,
            DPB1TCEStatus.UNEXPECTED_ERROR
        ], f"API v3.0 returned error status: {result}")

    def test_invalid_alleles_handling(self):
        """Test that API properly rejects invalid alleles with HTTP 422"""
        result = query_dpb1_tce(
            patient_dpb1="99:99",  # Invalid allele
            patient_dpb2="01:01",
            donor_dpb1="01:01",
            donor_dpb2="01:01",
            version="2.1",
            timeout=30
        )
        # Returns HTTP 422 which currently becomes REQUEST_ERROR
        self.assertEqual(
            result,
            DPB1TCEStatus.REQUEST_ERROR,
            "API should reject invalid alleles with HTTP 422"
        )


if __name__ == '__main__':
    # Allow running locally with:
    # RUN_INTEGRATION_TESTS=1 python -m pytest tests/test_integration.py
    os.environ['RUN_INTEGRATION_TESTS'] = '1'
    unittest.main()
