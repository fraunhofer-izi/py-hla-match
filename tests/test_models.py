import unittest
from py_hla_match.models import HLAPair, Individual, Patient, Donor
from py_hla_match.hla import HLA
from py_hla_match.exceptions import InvalidLocusComparisonError


class TestHLAPair(unittest.TestCase):

    def test_valid_hla_pair(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        self.assertEqual(pair.hla1, hla1)
        self.assertEqual(pair.hla2, hla2)
        self.assertEqual(pair.locus, "A")

    def test_hla_pair_with_one_none(self):
        hla1 = HLA("B*07:02")
        pair = HLAPair(hla1, None)
        self.assertEqual(pair.hla1, hla1)
        self.assertIsNone(pair.hla2)
        self.assertEqual(pair.locus, "B")

    def test_hla_pair_with_both_none(self):
        pair = HLAPair(None, None)
        self.assertIsNone(pair.hla1)
        self.assertIsNone(pair.hla2)
        self.assertIsNone(pair.locus)

    def test_hla_pair_with_invalid_type_raises_error(self):
        with self.assertRaises(TypeError):
            HLAPair("A*02:01", None)

    def test_hla_pair_with_mismatched_loci_raises_error(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("B*07:02")
        with self.assertRaises(InvalidLocusComparisonError):
            HLAPair(hla1, hla2)

    def test_hla_pair_with_drbx_handling(self):
        hla1 = HLA("DRB3*01:01")
        hla2 = HLA("DRB4*01:03")
        pair = HLAPair(hla1, hla2)
        self.assertEqual(pair.locus, "DRBX")

    def test_has_any_data_with_both_alleles(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        self.assertTrue(pair.has_any_data())

    def test_has_any_data_with_one_allele(self):
        hla1 = HLA("A*02:01")
        pair = HLAPair(hla1, None)
        self.assertTrue(pair.has_any_data())

    def test_has_any_data_with_no_alleles(self):
        pair = HLAPair(None, None)
        self.assertFalse(pair.has_any_data())

    def test_has_locus_info(self):
        hla1 = HLA("C*07:02")
        pair = HLAPair(hla1, None)
        self.assertTrue(pair.has_locus_info())

    def test_has_valid_hla_with_two_field_resolution(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        self.assertTrue(pair.has_valid_hla())

    def test_has_valid_hla_with_one_field_resolution(self):
        hla1 = HLA("A*02")
        pair = HLAPair(hla1, None)
        self.assertTrue(pair.has_valid_hla())

    def test_has_high_resolution_hla(self):
        hla1 = HLA("A*02:01")
        pair = HLAPair(hla1, None)
        self.assertTrue(pair.has_high_resolution_hla())

    def test_has_hla_pair_true(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        self.assertTrue(pair.has_hla_pair())

    def test_has_hla_pair_false(self):
        hla1 = HLA("A*02:01")
        pair = HLAPair(hla1, None)
        self.assertFalse(pair.has_hla_pair())

    def test_has_valid_hla_pair_both_valid(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        self.assertTrue(pair.has_valid_hla_pair())

    def test_has_valid_hla_pair_one_invalid(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*")  # Invalid - no allele group
        pair = HLAPair(hla1, hla2)
        self.assertFalse(pair.has_valid_hla_pair())

    def test_has_high_resolution_hla_pair_both_high_res(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        self.assertTrue(pair.has_high_resolution_hla_pair())

    def test_has_high_resolution_hla_pair_one_low_res(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03")  # Low resolution
        pair = HLAPair(hla1, hla2)
        self.assertFalse(pair.has_high_resolution_hla_pair())

    def test_str_representation(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        expected = "HLAPair(locus=A, hla1=A*02:01, hla2=A*03:01)"
        self.assertEqual(str(pair), expected)

    def test_str_representation_with_none_locus(self):
        pair = HLAPair(None, None)
        expected = "HLAPair(locus=None, hla1=None, hla2=None)"
        self.assertEqual(str(pair), expected)

    def test_repr(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        self.assertEqual(repr(pair), str(pair))


class TestIndividual(unittest.TestCase):

    def test_valid_individual(self):
        pair1 = HLAPair(HLA("A*02:01"), HLA("A*03:01"))
        pair2 = HLAPair(HLA("B*07:02"), HLA("B*08:01"))
        individual = Individual([pair1, pair2])
        self.assertEqual(len(individual.hla_data), 2)
        self.assertEqual(individual.hla_data[0], pair1)
        self.assertEqual(individual.hla_data[1], pair2)

    def test_individual_with_duplicate_locus_raises_error(self):
        pair1 = HLAPair(HLA("A*02:01"), HLA("A*03:01"))
        pair2 = HLAPair(HLA("A*24:02"), HLA("A*25:01"))
        with self.assertRaises(ValueError) as context:
            Individual([pair1, pair2])
        self.assertIn("Duplicate loci found: ['A']", str(context.exception))

    def test_individual_with_empty_hla_data(self):
        individual = Individual([])
        self.assertEqual(len(individual.hla_data), 0)

    def test_individual_with_none_locus_pairs(self):
        pair1 = HLAPair(None, None)
        pair2 = HLAPair(HLA("B*07:02"), HLA("B*08:01"))
        individual = Individual([pair1, pair2])
        self.assertEqual(len(individual.hla_data), 2)

    def test_get_hla_summary_all_high_resolution(self):
        pair1 = HLAPair(HLA("A*02:01"), HLA("A*03:01"))
        pair2 = HLAPair(HLA("B*07:02"), HLA("B*08:01"))
        individual = Individual([pair1, pair2])
        summary = individual.get_hla_summary()
        self.assertEqual(summary["total_alles"], 2)
        self.assertEqual(summary["parsed_alleles"], 2)
        self.assertEqual(summary["valid_alleles"], 2)
        self.assertEqual(summary["high_resolution_alleles"], 2)

    def test_get_hla_summary_mixed_resolution(self):
        pair1 = HLAPair(HLA("A*02:01"), HLA("A*03"))
        pair2 = HLAPair(HLA("B*07"), None)
        pair3 = HLAPair(None, None)
        individual = Individual([pair1, pair2, pair3])
        summary = individual.get_hla_summary()
        self.assertEqual(summary["total_alles"], 3)
        self.assertEqual(summary["parsed_alleles"], 2)
        self.assertEqual(summary["valid_alleles"], 1)
        self.assertEqual(summary["high_resolution_alleles"], 0)

    def test_get_hla_summary_empty_individual(self):
        individual = Individual([])
        summary = individual.get_hla_summary()
        self.assertEqual(summary["total_alles"], 0)
        self.assertEqual(summary["parsed_alleles"], 0)
        self.assertEqual(summary["valid_alleles"], 0)
        self.assertEqual(summary["high_resolution_alleles"], 0)


class TestPatientAndDonor(unittest.TestCase):

    def test_patient_creation(self):
        pair = HLAPair(HLA("A*02:01"), HLA("A*03:01"))
        patient = Patient([pair])
        self.assertIsInstance(patient, Patient)
        self.assertIsInstance(patient, Individual)
        self.assertEqual(len(patient.hla_data), 1)

    def test_donor_creation(self):
        pair = HLAPair(HLA("B*07:02"), HLA("B*08:01"))
        donor = Donor([pair])
        self.assertIsInstance(donor, Donor)
        self.assertIsInstance(donor, Individual)
        self.assertEqual(len(donor.hla_data), 1)

    def test_patient_inherits_individual_methods(self):
        pair1 = HLAPair(HLA("A*02:01"), HLA("A*03:01"))
        pair2 = HLAPair(HLA("B*07:02"), HLA("B*08:01"))
        patient = Patient([pair1, pair2])
        summary = patient.get_hla_summary()
        self.assertEqual(summary["total_alles"], 2)
        self.assertEqual(summary["high_resolution_alleles"], 2)

    def test_donor_inherits_individual_methods(self):
        pair1 = HLAPair(HLA("C*07:02"), HLA("C*08:01"))
        donor = Donor([pair1])
        summary = donor.get_hla_summary()
        self.assertEqual(summary["total_alles"], 1)
        self.assertEqual(summary["high_resolution_alleles"], 1)


if __name__ == "__main__":
    unittest.main()
