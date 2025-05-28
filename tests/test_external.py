import unittest
from unittest.mock import patch, Mock
import requests

from py_hla_match.external import query_dpb1_tce, DPB1TCEStatus


class TestDPB1TCE(unittest.TestCase):

    def create_mock_response(self, tce_prediction, version="2.1"):
        response_keys = {
            "2.0": "HLA-DPB1_TCE_report_V2.0",
            "2.1": "HLA-DPB1_TCE_report_V2.1",
            "3.0": "HLA-DPB1_TCE_report_V3.0"
        }

        return {
            response_keys[version]: {
                'donors': [{
                    'result': {
                        'tce_prediction': tce_prediction
                    }
                }]
            }
        }

    @patch('requests.get')
    def test_permissive_match(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Permissive")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "01:01", "02:01")
        self.assertEqual(result, DPB1TCEStatus.PERMISSIVE)

    @patch('requests.get')
    def test_non_permissive_gvh(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Non-Permissive GvH")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("10:01", "10:01", "14:01", "14:01")
        self.assertEqual(result, DPB1TCEStatus.NON_PERMISSIVE_GVH)

    @patch('requests.get')
    def test_non_permissive_hvg(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Non-Permissive HvG")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("14:01", "14:01", "10:01", "10:01")
        self.assertEqual(result, DPB1TCEStatus.NON_PERMISSIVE_HVG)

    @patch('requests.get')
    def test_invalid_alleles(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = self.create_mock_response(
            "Not possible as typing contains non-existent allele"
        )
        mock_get.return_value = mock_response

        result = query_dpb1_tce("99:99", "10:01", "14:01", "20:01")
        self.assertEqual(result, DPB1TCEStatus.INVALID_ALLELES)

    @patch('requests.get')
    def test_api_version_2_0(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Permissive", "2.0")
        mock_get.return_value = mock_response

        result = query_dpb1_tce(
            "01:01", "02:01", "03:01", "04:01", version="2.0"
        )
        self.assertEqual(result, DPB1TCEStatus.PERMISSIVE)
        # verify endpoint was called correctly
        mock_get.assert_called_once()
        call_args = mock_get.call_args
        self.assertIn("dpb1_tce_v2", call_args[0][0])

    @patch('requests.get')
    def test_api_version_2_1(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Permissive", "2.1")
        mock_get.return_value = mock_response

        result = query_dpb1_tce(
            "01:01", "02:01", "03:01", "04:01", version="2.1"
        )
        self.assertEqual(result, DPB1TCEStatus.PERMISSIVE)
        # verify endpoint was called correctly
        mock_get.assert_called_once()
        call_args = mock_get.call_args
        self.assertIn("dpb1_tce_v21", call_args[0][0])

    @patch('requests.get')
    def test_timeout_error(self, mock_get):
        mock_get.side_effect = requests.Timeout("Connection timed out")

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01", timeout=5)
        self.assertEqual(result, DPB1TCEStatus.TIMEOUT_ERROR)

    @patch('requests.get')
    def test_request_exception(self, mock_get):
        mock_get.side_effect = requests.RequestException("Network error")

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.REQUEST_ERROR)

    @patch('requests.get')
    def test_invalid_json_response(self, mock_get):
        mock_response = Mock()
        mock_response.json.side_effect = ValueError("Invalid JSON")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.VALUE_ERROR)

    @patch('requests.get')
    def test_missing_key_in_response(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = {"wrong_key": "wrong_value"}
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.API_ERROR)

    @patch('requests.get')
    def test_empty_donors_list(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = {
            "HLA-DPB1_TCE_report_V3.0": {
                'patient': {
                    'donors': []  # empty donors
                }
            }
        }
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.API_ERROR)

    @patch('requests.get')
    def test_missing_tce_prediction(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = {
            "HLA-DPB1_TCE_report_V3.0": {
                'patient': {
                    'donors': [{
                        'result': {}  # no tce_prediction
                    }]
                }
            }
        }
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.API_ERROR)

    @patch('requests.get')
    def test_unknown_tce_prediction(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Unknown Status")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.API_ERROR)

    def test_invalid_api_version(self):
        result = query_dpb1_tce(
            "01:01", "02:01", "03:01", "04:01", version="4.0"
        )
        self.assertEqual(result, DPB1TCEStatus.CONFIGURATION_ERROR)

    @patch('requests.get')
    def test_unexpected_exception(self, mock_get):
        mock_get.side_effect = Exception("Unexpected error")

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.UNEXPECTED_ERROR)

    @patch('requests.get')
    def test_query_parameters(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Permissive")
        mock_get.return_value = mock_response

        query_dpb1_tce("01:01", "02:01", "03:01", "04:01")

        # verify query parameters
        mock_get.assert_called_once()
        call_kwargs = mock_get.call_args[1]
        params = call_kwargs['params']
        self.assertEqual(params['pid'], 1)
        self.assertEqual(params['patdpb1'], '01:01')
        self.assertEqual(params['patdpb2'], '02:01')
        self.assertEqual(params['did'], 2)
        self.assertEqual(params['dondpb1'], '03:01')
        self.assertEqual(params['dondpb2'], '04:01')
        self.assertEqual(call_kwargs['timeout'], 10)

    @patch('requests.get')
    def test_custom_timeout(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Permissive")
        mock_get.return_value = mock_response

        query_dpb1_tce("01:01", "02:01", "03:01", "04:01", timeout=30)

        # verify timeout
        call_kwargs = mock_get.call_args[1]
        self.assertEqual(call_kwargs['timeout'], 30)

    @patch('requests.get')
    def test_http_error_status(self, mock_get):
        mock_response = Mock()
        mock_response.raise_for_status.side_effect = \
            requests.HTTPError("404 Not Found")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "03:01", "04:01")
        self.assertEqual(result, DPB1TCEStatus.REQUEST_ERROR)

    @patch('requests.get')
    def test_ard_matched_maps_to_permissive(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("ARD Matched")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "01:01", "02:01")
        self.assertEqual(result, DPB1TCEStatus.PERMISSIVE)

    @patch('requests.get')
    def test_permissive_core_maps_to_permissive(self, mock_get):
        mock_response = Mock()
        mock_response.json.return_value = \
            self.create_mock_response("Permissive (Core)")
        mock_get.return_value = mock_response

        result = query_dpb1_tce("01:01", "02:01", "01:01", "02:01")
        self.assertEqual(result, DPB1TCEStatus.PERMISSIVE)


if __name__ == '__main__':
    unittest.main()
