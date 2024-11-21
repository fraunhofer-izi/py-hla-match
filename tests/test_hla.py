import unittest
from py_hla_match.exceptions import MalformedHLAStringError
from py_hla_match import HLA

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

    def test_repr(self):
        hla = HLA("HLA-A*32:11Q")
        expected_repr = (
            "HLA(allele_string='HLA-A*32:11Q', locus='A', allele_group='32', "
            "allele='11', synonymous_variant=None, non_coding_variant=None, suffix='Q')"
        )
        self.assertEqual(repr(hla), expected_repr)
