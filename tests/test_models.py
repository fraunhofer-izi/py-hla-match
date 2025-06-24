import unittest
from py_hla_match.models import HLAPair, Individual, Patient, Donor
from py_hla_match.hla import HLA
from py_hla_match.exceptions import InvalidLocusComparisonError


class TestHLAPair(unittest.TestCase):
    """Tests for HLAPair."""

    def test_initialization_success(self):
        hla1 = HLA("A*01:01")
        hla2 = HLA("A*02:01")
        pair = HLAPair(hla1, hla2)
        self.assertEqual(pair.hla1, hla1)
        self.assertEqual(pair.hla2, hla2)
        self.assertEqual(pair.locus, "A")

    def test_initialization_raises_type_error_for_non_hla_object(self):
        hla1 = HLA("A*01:01")
        with self.assertRaisesRegex(TypeError, "must be an instance of HLA"):
            HLAPair(hla1, "not-an-hla-object")
        with self.assertRaisesRegex(TypeError, "must be an instance of HLA"):
            HLAPair("not-an-hla-object", hla1)

    def test_initialization_raises_error_for_mismatched_loci(self):
        hla1 = HLA("A*01:01")
        hla2 = HLA("B*07:02")
        with self.assertRaises(InvalidLocusComparisonError):
            HLAPair(hla1, hla2)

    def test_get_paired_resolution_is_minimum_of_pair(self):
        # 2-field and 4-field -> 2-field
        hla_2_field = HLA("A*01:01")
        hla_4_field = HLA("A*02:01:01:01")
        pair1 = HLAPair(hla_2_field, hla_4_field)
        self.assertEqual(pair1.get_paired_resolution(), 2)

        # 1-field and 3-field -> 1-field
        hla_1_field = HLA("B*07")
        hla_3_field = HLA("B*08:01:02")
        pair2 = HLAPair(hla_1_field, hla_3_field)
        self.assertEqual(pair2.get_paired_resolution(), 1)

    def test_equality_is_order_insensitive(self):
        hla1 = HLA("C*01:02")
        hla2 = HLA("C*02:03")
        pair1 = HLAPair(hla1, hla2)
        pair2 = HLAPair(hla2, hla1)
        self.assertEqual(pair1, pair2)

    def test_inequality_for_different_alleles(self):
        hla1 = HLA("A*01:01")
        hla2 = HLA("A*02:01")
        hla3 = HLA("A*03:01")
        pair1 = HLAPair(hla1, hla2)
        pair2 = HLAPair(hla1, hla3)
        self.assertNotEqual(pair1, pair2)

    def test_hashability_and_set_uniqueness(self):
        hla1 = HLA("DRB1*01:01")
        hla2 = HLA("DRB1*04:02")
        pair1 = HLAPair(hla1, hla2)
        pair2 = HLAPair(hla2, hla1)
        pair3 = HLAPair(hla1, hla1)

        # set should treat pair1 and pair2 as equal items
        hla_set = {pair1, pair2, pair3}
        self.assertEqual(len(hla_set), 2)
        self.assertIn(pair1, hla_set)
        self.assertIn(pair3, hla_set)

    def test_str_and_repr_representation(self):
        hla1 = HLA("A*02:01")
        hla2 = HLA("A*03:01")
        pair = HLAPair(hla1, hla2)
        expected_str = "HLAPair(locus=A, hla1=A*02:01, hla2=A*03:01)"
        self.assertEqual(str(pair), expected_str)
        self.assertEqual(repr(pair), expected_str)


class TestDRB345Pairs(unittest.TestCase):
    """Edge-cases specific to the DRB3/4/5/X region."""
    def test_pair_drb3_drb3(self):
        pair = HLAPair(HLA("DRB3*01:01"), HLA("DRB3*02:02"))
        self.assertEqual(pair.locus, "DRB345")
        self.assertEqual(pair.get_paired_resolution(), 2)

    def test_pair_drb3_drb4(self):
        # different genes, same canonical locus -> still a legal pair
        pair = HLAPair(HLA("DRB3*01:01"), HLA("DRB4*01:01"))
        self.assertEqual(pair.locus, "DRB345")

    def test_pair_drb3_missing(self):
        pair = HLAPair(HLA("DRB3*01:01"), HLA("DRBX*NE"))
        self.assertEqual(pair.locus, "DRB345")
        # one ‘NE’ allele => min-resolution 0
        self.assertEqual(pair.get_paired_resolution(), 0)

    def test_pair_drb3_with_drb1_raises(self):
        with self.assertRaises(InvalidLocusComparisonError):
            HLAPair(HLA("DRB3*01:01"), HLA("DRB1*01:01"))

    def test_duplicate_drb345_pairs_raise(self):
        pair1 = HLAPair(HLA("DRB3*01:01"), HLA("DRBX*NE"))
        pair2 = HLAPair(HLA("DRB4*01:01"), HLA("DRB5*01:01"))
        with self.assertRaisesRegex(ValueError, "Multiple .*'DRB345'.* found"):
            Individual([pair1, pair2])

    def test_individual_single_drb345_pair(self):
        pair = HLAPair(HLA("DRB4*01:01"), HLA("DRB3*02:02"))
        ind = Individual([pair])
        self.assertEqual(len(ind.hla_data), 1)
        self.assertEqual(ind.hla_data[0].locus, "DRB345")


class TestIndividual(unittest.TestCase):
    """Tests for Individual."""

    def test_initialization_success(self):
        pair1 = HLAPair(HLA("A*01:01"), HLA("A*02:01"))
        pair2 = HLAPair(HLA("B*07:02"), HLA("B*08:01"))
        individual = Individual([pair1, pair2])
        self.assertEqual(len(individual.hla_data), 2)
        self.assertIn(pair1, individual.hla_data)

    def test_initialization_with_empty_list(self):
        individual = Individual([])
        self.assertEqual(len(individual.hla_data), 0)

    def test_initialization_raises_type_error_for_invalid_item_in_list(self):
        pair1 = HLAPair(HLA("A*01:01"), HLA("A*02:01"))
        invalid_item = "not-a-pair"
        with self.assertRaisesRegex(
            TypeError, "must contain only HLAPair objects"
        ):
            Individual([pair1, invalid_item])

    def test_initialization_raises_value_error_for_duplicate_loci(self):
        pair1 = HLAPair(HLA("A*01:01"), HLA("A*02:01"))
        pair2 = HLAPair(HLA("A*24:02"), HLA("A*25:01"))
        with self.assertRaisesRegex(ValueError, "Multiple .*'A'.* found"):
            Individual([pair1, pair2])

    def test_get_hla_summary_with_data(self):
        pair_a = HLAPair(HLA("A*01:01"), HLA("A*02:01:01"))  # Resolution -> 2
        pair_b = HLAPair(HLA("B*07:02"), HLA("B*08:01"))  # Resolution -> 2
        pair_c = HLAPair(HLA("C*01"), HLA("C*01:02:01:01"))  # Resolution -> 1
        individual = Individual([pair_a, pair_b, pair_c])

        summary = individual.get_hla_summary()

        expected_summary = {
            "total_loci_typed": 3,
            "resolution_summary": {2: 2, 1: 1}  # 2 x 2-field, 1 x 1-field
        }
        self.assertEqual(summary, expected_summary)

    def test_get_hla_summary_on_empty_individual(self):
        individual = Individual([])
        summary = individual.get_hla_summary()
        expected_summary = {
            "total_loci_typed": 0,
            "resolution_summary": {}
        }
        self.assertEqual(summary, expected_summary)


class TestPatientAndDonor(unittest.TestCase):
    """Test Patient and Donor classes inherit correctly from Individual."""

    def test_patient_and_donor_creation(self):
        pair = HLAPair(HLA("A*01:01"), HLA("A*02:01"))
        patient = Patient([pair])
        donor = Donor([pair])

        self.assertIsInstance(patient, Patient)
        self.assertIsInstance(patient, Individual)
        self.assertIsInstance(donor, Donor)
        self.assertIsInstance(donor, Individual)

    def test_patient_and_donor_inherit_individual_behavior(self):
        # sanity check inheritance
        pair1 = HLAPair(HLA("A*01:01"), HLA("A*02:01"))
        pair2 = HLAPair(HLA("A*24:02"), HLA("A*25:01"))
        with self.assertRaises(ValueError):
            Patient([pair1, pair2])

        # get_hla_summary inheritance
        donor = Donor([pair1])
        summary = donor.get_hla_summary()
        self.assertEqual(summary["total_loci_typed"], 1)
        self.assertEqual(summary["resolution_summary"], {2: 1})


if __name__ == "__main__":
    unittest.main()
