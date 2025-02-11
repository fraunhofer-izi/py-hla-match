import unittest
from py_hla_match.exceptions import MalformedHLAStringError
from py_hla_match.hla import HLA


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

    def test_invalid_suffix(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("HLA-A*32:11X")

    def test_invalid_suffix_numeric(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("HLA-A*24:02:01:01:01")

    def test_invalid_suffix_special_character(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("HLA-B*44:02@")

    def test_invalid_suffix_lowercase(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("HLA-A*24:02l")

    def test_invalid_hla_with_incorrect_field(self):
        with self.assertRaises(MalformedHLAStringError):
            HLA("HLA-A*02:01X:01")

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
