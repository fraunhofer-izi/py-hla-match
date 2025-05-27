import unittest

from py_hla_match.exceptions import MalformedHLAStringError
from py_hla_match.hla import HLA
import pyard


class TestHLA(unittest.TestCase):

    def test_valid_hla_with_suffix(self):
        hla = HLA("HLA-A*32:11Q")
        self.assertEqual(hla.locus, "A")
        self.assertEqual(hla.allele_group, "32")
        self.assertEqual(hla.allele, "11")
        self.assertIsNone(hla.synonymous_variant)
        self.assertIsNone(hla.non_coding_variant)
        self.assertEqual(hla.suffix, "Q")

    def test_valid_hla_without_suffix(self):
        hla = HLA("B*07:02")
        self.assertEqual(hla.locus, "B")
        self.assertEqual(hla.allele_group, "07")
        self.assertEqual(hla.allele, "02")
        self.assertIsNone(hla.synonymous_variant)
        self.assertIsNone(hla.non_coding_variant)
        self.assertIsNone(hla.suffix)

    def test_valid_hla_with_prefix_removed(self):
        hla = HLA("HLA-DRB1*15:01")
        self.assertEqual(hla.locus, "DRB1")
        self.assertEqual(hla.allele_group, "15")
        self.assertEqual(hla.allele, "01")
        self.assertIsNone(hla.synonymous_variant)
        self.assertIsNone(hla.non_coding_variant)
        self.assertIsNone(hla.suffix)

    def test_valid_hla_with_synonymous_variant_only(self):
        hla = HLA("HLA-DQB1*03:01:02")
        self.assertEqual(hla.locus, "DQB1")
        self.assertEqual(hla.allele_group, "03")
        self.assertEqual(hla.allele, "01")
        self.assertEqual(hla.synonymous_variant, "02")
        self.assertIsNone(hla.non_coding_variant)
        self.assertIsNone(hla.suffix)

    def test_valid_hla_with_non_coding_variant(self):
        hla = HLA("HLA-DRB1*13:01:01:02")
        self.assertEqual(hla.locus, "DRB1")
        self.assertEqual(hla.allele_group, "13")
        self.assertEqual(hla.allele, "01")
        self.assertEqual(hla.synonymous_variant, "01")
        self.assertEqual(hla.non_coding_variant, "02")
        self.assertIsNone(hla.suffix)

    def test_valid_hla_with_null_suffix(self):
        hla = HLA("HLA-A*24:09N")
        self.assertEqual(hla.locus, "A")
        self.assertEqual(hla.allele_group, "24")
        self.assertEqual(hla.allele, "09")
        self.assertIsNone(hla.synonymous_variant)
        self.assertIsNone(hla.non_coding_variant)
        self.assertEqual(hla.suffix, "N")

    def test_valid_hla_with_low_expression_suffix(self):
        hla = HLA("HLA-A*30:14L")
        self.assertEqual(hla.locus, "A")
        self.assertEqual(hla.allele_group, "30")
        self.assertEqual(hla.allele, "14")
        self.assertIsNone(hla.synonymous_variant)
        self.assertIsNone(hla.non_coding_variant)
        self.assertEqual(hla.suffix, "L")

    def test_valid_hla_with_secreted_suffix(self):
        hla = HLA("HLA-B*44:02:01:02S")
        self.assertEqual(hla.locus, "B")
        self.assertEqual(hla.allele_group, "44")
        self.assertEqual(hla.allele, "02")
        self.assertEqual(hla.synonymous_variant, "01")
        self.assertEqual(hla.non_coding_variant, "02")
        self.assertEqual(hla.suffix, "S")

    def test_valid_hla_without_prefix(self):
        hla = HLA("A*01:01:01")
        self.assertEqual(hla.locus, "A")
        self.assertEqual(hla.allele_group, "01")
        self.assertEqual(hla.allele, "01")
        self.assertEqual(hla.synonymous_variant, "01")
        self.assertIsNone(hla.non_coding_variant)
        self.assertIsNone(hla.suffix)
        self.assertIsNone(hla.group_code)

    def test_valid_hla_with_group_code_G(self):
        hla = HLA("DQB1*06:02:01G")
        self.assertEqual(hla.locus, "DQB1")
        self.assertEqual(hla.allele_group, "06")
        self.assertEqual(hla.allele, "02")
        self.assertEqual(hla.synonymous_variant, "01")
        self.assertIsNone(hla.non_coding_variant)
        self.assertIsNone(hla.suffix)
        self.assertEqual(hla.group_code, "G")

    def test_valid_hla_with_group_code_P(self):
        hla = HLA("DQA1*05:88P")
        self.assertEqual(hla.locus, "DQA1")
        self.assertEqual(hla.allele_group, "05")
        self.assertEqual(hla.allele, "88")
        self.assertIsNone(hla.synonymous_variant)
        self.assertIsNone(hla.non_coding_variant)
        self.assertIsNone(hla.suffix)
        self.assertEqual(hla.group_code, "P")

    def test_valid_hla_with_three_digit_fields(self):
        hla = HLA("A*68:114:02")
        self.assertEqual(hla.locus, "A")
        self.assertEqual(hla.allele_group, "68")
        self.assertEqual(hla.allele, "114")
        self.assertEqual(hla.synonymous_variant, "02")
        self.assertIsNone(hla.non_coding_variant)
        self.assertIsNone(hla.suffix)
        self.assertIsNone(hla.group_code)

    def test_valid_hla_with_full_fields(self):
        hla = HLA("C*06:127:01:02")
        self.assertEqual(hla.locus, "C")
        self.assertEqual(hla.allele_group, "06")
        self.assertEqual(hla.allele, "127")
        self.assertEqual(hla.synonymous_variant, "01")
        self.assertEqual(hla.non_coding_variant, "02")
        self.assertIsNone(hla.suffix)
        self.assertIsNone(hla.group_code)

    def test_invalid_locus(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("HLA-Ad*32:11")

    def test_missing_locus(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("*24:02:01:01")

    def test_missing_asterix(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("HLA-B44:02")

    def test_invalid_allele(self):
        with self.assertRaises(pyard.exceptions.InvalidAlleleError):
            HLA("HLA-DPB1*07:32")

    def test_incomplete_allele_logs_warning(self):
        with self.assertLogs('py_hla_match.hla', level='WARNING') as captured:
            HLA("HLA-A*32")

        self.assertTrue(
            any(
                "is not a specific allele." in message for message
                in captured.output
            ),
            "Expected the log message about 'not a specific allele.'"
        )
        self.assertTrue(
            any(
                "Validity of '32' is not checked." in message for message
                in captured.output
            ),
            "Expected the log message about 'Validity of ... is not checked.'"
        )

    def test_logger_warnings_for_unparseable_content(self):
        with self.assertLogs('py_hla_match.hla', level='WARNING') as captured:
            HLA("HLA-C**invalid")

        # Now check that the expected warnings appeared
        self.assertTrue(
            any(
                "did match HLA Nomenclature" in msg for msg
                in captured.output
            ),
            "Expected 'did match HLA Nomenclature' warning in logs"
        )
        self.assertTrue(
            any("has unparseable content" in msg for msg in captured.output),
            "Expected 'has unparseable content' warning in logs"
        )

    def test_repr(self):
        hla = HLA("HLA-A*32:11Q")
        expected_repr = (
            "HLA(allele_string='HLA-A*32:11Q', locus='A', allele_group='32', "
            "allele='11', synonymous_variant=None, non_coding_variant=None, "
            "suffix='Q', group_code=None, "
            "ard_redux_allele_string='HLA-A*32:11Q', "
            "ard_redux_allele_group='32', ard_redux_allele='11')"
        )
        self.assertEqual(repr(hla), expected_repr)


if __name__ == "__main__":
    unittest.main()
