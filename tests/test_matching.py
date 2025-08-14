import unittest

from py_hla_match.hla import HLA
from py_hla_match.matching import (
    allele_match,
    _get_correct_allele_pairing,
    allele_pair_match,
    AlleleMatchLevel,
    MatchResult,
    multi_locus_match
)
from py_hla_match.exceptions import (
    InvalidLocusComparisonError, MalformedHLAStringError
)
from pyard.exceptions import InvalidAlleleError
from py_hla_match.models import HLAPair, Individual


class TestAlleleMatch(unittest.TestCase):
    def test_invalid_hla_object(self):
        # Test with hla1 None
        allele1 = None
        allele2 = HLA("A*01:01")
        with self.assertRaises(TypeError) as context:
            allele_match(allele1, allele2)
        self.assertIn(
            "hla1 must be an instance of HLA", str(context.exception)
        )

        # Test with hla2 as a string
        allele1 = HLA("A*01:01")
        allele2 = "Not an HLA object"
        with self.assertRaises(TypeError) as context:
            allele_match(allele1, allele2)
        self.assertIn(
            "hla2 must be an instance of HLA", str(context.exception)
        )

        # Test with both hla1 and hla2 invalid (integers)
        allele1 = 123
        allele2 = 456
        with self.assertRaises(TypeError) as context:
            allele_match(allele1, allele2)
        self.assertIn(
            "hla1 must be an instance of HLA", str(context.exception)
        )

    def test_valid_drb34_locus_mismatch(self):
        """
        Test Case: Locus mismatch
        Allele 1: DRB3*02:02:01
        Allele 2: DRB4*01:03:01
        Expected Match Level: LOCUS_MISMATCH (0)
        """
        allele1 = HLA("DRB3*02:02:01")
        allele2 = HLA("DRB4*01:03:01")
        expected_match_level = AlleleMatchLevel.LOCUS_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_drb35_locus_mismatch(self):
        """
        Test Case: Locus mismatch
        Allele 1: DRB3*01
        Allele 2: DRB5*01
        Expected Match Level: LOCUS_MISMATCH (0)
        """
        allele1 = HLA("DRB3*01")
        allele2 = HLA("DRB5*01")
        expected_match_level = AlleleMatchLevel.LOCUS_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_drb3x_locus_mismatch(self):
        """
        Test Case: Locus mismatch
        Allele 1: DRB3*01
        Allele 2: DRBX*NE
        Expected Match Level: LOCUS_MISMATCH (0)
        """
        allele1 = HLA("DRB3*01")
        allele2 = HLA("DRBX*NE")
        expected_match_level = AlleleMatchLevel.LOCUS_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_invalid_locus_mismatch(self):
        """
        Test Case: Invalid locus mismatch
        InvalidLocusComparisonError
        Allele 1: A*01:01:01:01
        Allele 2: B*07:02:01:01
        Expected: InvalidLocusComparisonError
        """
        allele1 = HLA("A*01:01:01:01")
        allele2 = HLA("B*07:02:01:01")
        with self.assertRaises(InvalidLocusComparisonError) as context:
            allele_match(allele1, allele2)
        expected_message = (
            f"Invalid locus comparison between '{allele1.locus}' and "
            f"'{allele2.locus}'."
        )
        self.assertEqual(str(context.exception), expected_message)

    def test_valid_allele_group_mismatch(self):
        """
        Test Case: Allele group mismatch
        Allele 1: DPB1*02:01:02
        Allele 2: DPB1*04:02:01
        Expected Match Level: ALLELE_GROUP_MISMATCH (1)
        """
        allele1 = HLA("DPB1*02:01:02")
        allele2 = HLA("DPB1*04:02:01")
        expected_match_level = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_allele_mismatch_with_P_group(self):
        """
        Test Case: Allele mismatch with 'P' suffix
        Allele 1: DPB1*04:01P
        Allele 2: DPB1*04:02P
        Expected Match Level: ALLELE_MISMATCH (2)
        """
        allele1 = HLA("DPB1*04:01P")
        allele2 = HLA("DPB1*04:02P")
        expected_match_level = AlleleMatchLevel.ALLELE_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_ard_match_without_group_code(self):
        """
        Test Case: ARD match without group code
        Allele 1: A*02:01
        Allele 2: A*02:1193
        Expected Match Level: ARD_MATCH (3)
        """
        allele1 = HLA("A*02:01")
        allele2 = HLA("A*02:1193")
        expected_match_level = AlleleMatchLevel.ARD_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_ard_match_with_G_group(self):
        """
        Test Case: ARD match with G group
        Allele 1: C*07:02:01G
        Allele 2: C*07:1058
        Expected Match Level: ARD_MATCH (3)
        """
        allele1 = HLA("C*07:02:01G")
        allele2 = HLA("C*07:1058")
        expected_match_level = AlleleMatchLevel.ARD_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_not_applicable_with_Q_suffix(self):
        """
        Test Case: ARD match with Q suffix
        Allele 1: A*01:436Q
        Allele 2: A*01:01:70
        Expected Match Level: NOT_APPLICABLE
        """
        allele1 = HLA("A*01:436Q")
        allele2 = HLA("A*01:01:70")
        expected_match_level = AlleleMatchLevel.NOT_APPLICABLE
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_allele_mismatch_with_L_suffix(self):
        """
        Test Case: ARD match with L suffix
        Allele 1: B*38:68L
        Allele 2: B*38:01P
        Expected Match Level: ALLELE_MISMATCH (2)
        """
        allele1 = HLA("B*38:68L")
        allele2 = HLA("B*38:01P")
        expected_match_level = AlleleMatchLevel.ALLELE_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_allele_mismatch_with_N_suffix(self):
        """
        Test Case: ARD match with N suffix
        Allele 1: C*03:693
        Allele 2: C*03:20N
        Expected Match Level: ALLELE_MISMATCH (2)
        """
        # TODO: external validation of correct logic required
        # currently resolved to ARD_MATCH, solely relying on py-ard
        allele1 = HLA("C*03:693")
        allele2 = HLA("C*03:20N")
        expected_match_level = AlleleMatchLevel.ALLELE_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_synonymous_variant_match(self):
        """
        Test Case: Synonymous variant match with missing 4-field
        Allele 1: DPA1*01:03:01
        Allele 2: DPA1*01:03:01
        Expected Match Level: SYNONYMOUS_VARIANT_MATCH (4)
        """
        allele1 = HLA("DPA1*01:03:01")
        allele2 = HLA("DPA1*01:03:01")
        expected_match_level = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_non_encoding_variant_match(self):
        """
        Test Case: Non encoding variant match
        Allele 1: A*01:01:01:46
        Allele 2: A*01:01:01:46
        Expected Match Level: NON_ENCODING_VARIANT_MATCH (5)
        """
        allele1 = HLA("A*01:01:01:46")
        allele2 = HLA("A*01:01:01:46")
        expected_match_level = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_not_applicable_due_to_nan(self):
        """
        Test Case: missing allele data *NE
        Allele 1: A*NE
        Allele 2: A*01:01
        Expected Match Level: NOT_APPLICABLE
        """
        allele1 = HLA("A*NE")
        allele2 = HLA("A*01:01")
        expected_match_level = AlleleMatchLevel.NOT_APPLICABLE
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_one_field_vs_two_field_same_group_not_applicable(self):
        """
        Test Case: 1-field vs 2-field with same group
        Allele 1: B*07
        Allele 2: B*07:05
        Expected Match Level: NOT_APPLICABLE
        """
        allele1 = HLA("B*07")
        allele2 = HLA("B*07:05")
        expected_match_level = AlleleMatchLevel.NOT_APPLICABLE
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_allele_group_mismatch_with_one_field(self):
        """
        Test Case: 1-field vs 1-field with different groups
        Allele 1: C*01
        Allele 2: C*02
        Expected Match Level: ALLELE_GROUP_MISMATCH
        """
        allele1 = HLA("C*01")
        allele2 = HLA("C*02")
        expected_match_level = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_not_applicable_with_suffix_Q(self):
        """
        Test Case: Suffix does match beyond ARD_MATCH
        Allele 1: A*24:473Q
        Allele 2: A*02:99
        Expected Match Level: ARD_MATCH
        """
        allele1 = HLA("A*24:473Q")
        allele2 = HLA("A*24:02P")
        expected_match_level = AlleleMatchLevel.NOT_APPLICABLE
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_ard_match_with_different_group_codes(self):
        """
        Test Case: ARD match with G vs P group-codes
        Allele 1: A*01:02P
        Allele 2: A*01:02:01G
        Expected Match Level: ARD_MATCH
        """
        allele1 = HLA("A*01:02P")
        allele2 = HLA("A*01:02:01G")
        expected_match_level = AlleleMatchLevel.ARD_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_synonymous_vs_full_variant(self):
        """
        Test Case: 3-field vs 4-field with equal synonymous variant
        Allele 1: DQA1*01:01:02
        Allele 2: DQA1*01:01:02:07
        Expected Match Level: SYNONYMOUS_VARIANT_MATCH
        """
        allele1 = HLA("DQA1*01:01:02")
        allele2 = HLA("DQA1*01:01:02:07")
        expected_match_level = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)


class TestGetCorrectAllelePairing(unittest.TestCase):
    def test_tie_on_score_keeps_first_pairing(self):
        """
        Test Case: Tie on pairing score
        Patient Alleles: A*01:01:01:01, A*02:01
        Donor Alleles:   A*01:01:01:01, A*01:01:01:01
        Expected Pairing Levels: (allele_match(p1, d1), allele_match(p2, d2))
        Expected Best Score: sum(expected_pairing_levels)
        """
        p1, p2 = HLA("A*01:01:01:01"), HLA("A*02:01")
        d1, d2 = HLA("A*01:01:01:01"), HLA("A*01:01:01:01")
        patient = HLAPair(p1, p2)
        donor = HLAPair(d1, d2)

        best_score, chosen_pairing = \
            _get_correct_allele_pairing(patient, donor)

        expected_pairing = (
            allele_match(p1, d1),
            allele_match(p2, d2),
        )
        self.assertEqual(chosen_pairing, expected_pairing)
        self.assertEqual(best_score, sum(expected_pairing))

    def test_all_not_applicable_pairing_score(self):
        """
        Test Case: NOT_APPLICABLE twice
        Patient Alleles: B*NA, B*NA
        Donor Alleles:   B*NA, B*NA
        Expected Best Score: AlleleMatchLevel.NOT_APPLICABLE * 2
        Expected Pairing Levels:
            (AlleleMatchLevel.NOT_APPLICABLE, AlleleMatchLevel.NOT_APPLICABLE)
        """
        p1 = p2 = HLA("B*NA")
        d1 = d2 = HLA("B*NA")
        patient = HLAPair(p1, p2)
        donor = HLAPair(d1, d2)

        best_score, levels = _get_correct_allele_pairing(patient, donor)

        self.assertEqual(best_score, AlleleMatchLevel.NOT_APPLICABLE * 2)
        self.assertEqual(
            levels,
            (AlleleMatchLevel.NOT_APPLICABLE, AlleleMatchLevel.NOT_APPLICABLE),
        )

    def test_pairing_prefers_lower_negative_penalty(self):
        """
        Test Case: Two possible pairings – one sums to 0, the other to +1.
        """
        p = HLAPair(HLA("B*07"),     HLA("B*07:02"))
        d = HLAPair(HLA("B*07:02"),  HLA("B*07"))

        best_score, levels = _get_correct_allele_pairing(p, d)
        self.assertEqual(best_score, 1)
        self.assertEqual(
            levels,
            (AlleleMatchLevel.NOT_APPLICABLE, AlleleMatchLevel.ARD_MATCH)
        )

    def test_correct_pairing_with_ambiguous_alleles(self):
        """
        Test Case: pairing prefers NOT_APPLICABLE over any mismatch
        """
        p = HLAPair(HLA("B*07"),     HLA("B*01"))
        d = HLAPair(HLA("B*07"),  HLA("B*01"))

        best_score, levels = _get_correct_allele_pairing(p, d)
        self.assertEqual(best_score, 0)
        self.assertEqual(
            levels,
            (AlleleMatchLevel.NOT_APPLICABLE, AlleleMatchLevel.NOT_APPLICABLE)
        )


class TestAllelePairMatch(unittest.TestCase):
    def test_valid_match_without_swapping(self):
        """
        Test Case: Double SYNONYMOUS_VARIANT_MATCH (4)
        Patient: DRB1*15:01:01, DRB1*15:01:01
        Donor: DRB1*15:01:01, DRB1*15:01:01
        Expected Score: 8
        Expected Allele Match Levels: Double SYNONYMOUS_VARIANT_MATCH (4)
        """
        patient_allele1 = HLA("DRB1*15:01:01")
        patient_allele2 = HLA("DRB1*15:01:01")
        donor_allele1 = HLA("DRB1*15:01:01")
        donor_allele2 = HLA("DRB1*15:01:01")

        patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
        donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        result = allele_pair_match(patient, donor)

        expected_score = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH * 2
        expected_match_levels = (
            AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH,
            AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        )

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_match_levels)

    def test_valid_match_with_suffixes_G_and_P(self):
        """
        Test Case: Match with suffixes 'G' and 'P'
        Patient: DRB1*01:01:01G, DRB1*07:01:01G
        Donor: DRB1*01:01P, DRB1*07:01:01
        Expected Score: 6
        Expected Allele Match Levels: Double ARD_MATCH (3)
        """
        patient_allele1 = HLA("DRB1*01:01:01G")
        patient_allele2 = HLA("DRB1*07:01:01G")
        donor_allele1 = HLA("DRB1*01:01P")
        donor_allele2 = HLA("DRB1*07:01:01")

        patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
        donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        result = allele_pair_match(patient, donor)

        expected_score = AlleleMatchLevel.ARD_MATCH * 2
        expected_match_levels = (
            AlleleMatchLevel.ARD_MATCH,
            AlleleMatchLevel.ARD_MATCH
        )

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_match_levels)

    def test_valid_match_with_swapping(self):
        """
        Test Case: Match requiring swapping
        Patient: B*35:02:01, B*51:01P
        Donor: B*51:01:01, B*35:02:01
        Expected Score: 4 + 3
        Expected Allele Match Levels:
            - B*35:02:01 with B*35:02:01: SYNONYMOUS_VARIANT_MATCH
            - B*51:01P with B*51:01:01: ARD_MATCH
        """
        patient_allele1 = HLA("B*35:02:01")
        patient_allele2 = HLA("B*51:01P")
        donor_allele1 = HLA("B*51:01:01")
        donor_allele2 = HLA("B*35:02:01")

        patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
        donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        result = allele_pair_match(patient, donor)

        expected_match_levels = (
            AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH,
            AlleleMatchLevel.ARD_MATCH
        )
        expected_score = sum(expected_match_levels)

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(
            result.allele_match_levels,
            expected_match_levels
        )

    def test_valid_double_allele_mismatch(self):
        """
        Test Case: Double ALLELE_MISMATCH (2)
        Patient Alleles: DPB1*04:01:01, DPB1*04:01:01
        Donor Alleles: DPB1*04:02:01, DPB1*04:02:01
        Expected Score: 4
        Expected Allele Match Levels: Double ALLELE_MISMATCH (2)
        """
        patient_allele1 = HLA("DPB1*04:01:01")
        patient_allele2 = HLA("DPB1*04:01:01")
        donor_allele1 = HLA("DPB1*04:02:01")
        donor_allele2 = HLA("DPB1*04:02:01")

        patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
        donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        result = allele_pair_match(patient, donor)

        expected_score = AlleleMatchLevel.ALLELE_MISMATCH * 2
        expected_match_levels = (
            AlleleMatchLevel.ALLELE_MISMATCH,
            AlleleMatchLevel.ALLELE_MISMATCH
        )

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_match_levels)

    def test_valid_allele_group_mismatch_and_allele_mismatch_swapping(self):
        """
        Test Case: Swapping ALLELE_GROUP_MISMATCH, ALLELE_MISMATCH
        Patient Alleles: DPB1*01:01:01, DPB1*04:02:01
        Donor Alleles: DPB1*04:01:01, DPB1*02:01:02
        Expected Score: 3
        Expected Allele Match Levels: [ALLELE_MISMATCH, ALLELE_GROUP_MISMATCH]
        """
        patient_allele1 = HLA("DPB1*01:01:01")
        patient_allele2 = HLA("DPB1*04:02:01")
        donor_allele1 = HLA("DPB1*04:01:01")
        donor_allele2 = HLA("DPB1*02:01:02")

        patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
        donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        result = allele_pair_match(patient, donor)

        expected_match_levels = (
            AlleleMatchLevel.ALLELE_GROUP_MISMATCH,
            AlleleMatchLevel.ALLELE_MISMATCH
        )
        expected_score = sum(expected_match_levels)

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_match_levels)

    def test_valid_ard_match_and_allele_mismatch_swapping(self):
        """
        Test Case: ARD_MATCH and ALLELE_MISMATCH with swapping
        Patient Alleles: DPB1*03:01P, DPB1*04:01P
        Donor Alleles: 5
        Expected Allele Match Levels: [ARD_MATCH, ALLELE_MISMATCH]
        """
        patient_allele1 = HLA("DPB1*03:01P")
        patient_allele2 = HLA("DPB1*04:01P")
        donor_allele1 = HLA("DPB1*04:02P")
        donor_allele2 = HLA("DPB1*03:01P")

        patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
        donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        result = allele_pair_match(patient, donor)

        expected_match_levels = (
            AlleleMatchLevel.ARD_MATCH,
            AlleleMatchLevel.ALLELE_MISMATCH
        )
        expected_score = sum(expected_match_levels)

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_match_levels)

    def test_homozygous_patient_flag_true(self):
        """
        Test Case: Patient homozygous
        Patient Alleles: C*03:04:01, C*03:04:01
        Donor Alleles:   C*03:04:01, C*04:01:01
        Expected: is_homozygous_patient == True
        """
        allele1 = HLA("C*03:04:01")
        allele2 = HLA("C*03:04:01")
        allele3 = HLA("C*03:04:01")
        allele4 = HLA("C*04:01:01")
        patient = HLAPair(allele1, allele2)
        donor = HLAPair(allele3, allele4)

        result = allele_pair_match(patient, donor)
        self.assertTrue(result.is_homozygous_patient)

    def test_negative_score_when_not_applicable(self):
        """
        Test Case: Both comparisons NOT_APPLICABLE
        Patient Alleles: E*NE, E*NE
        Donor Alleles:   E*NE, E*NE
        Expected Pairing Score: AlleleMatchLevel.NOT_APPLICABLE * 2
        Expected Match Levels: (NOT_APPLICABLE, NOT_APPLICABLE)
        """
        na1 = HLA("E*NE")
        na2 = HLA("E*NE")
        na3 = HLA("E*NE")
        na4 = HLA("E*NE")
        patient = HLAPair(na1, na2)
        donor = HLAPair(na3, na4)

        result = allele_pair_match(patient, donor)

        expected_score = AlleleMatchLevel.NOT_APPLICABLE * 2
        expected_levels = (
            AlleleMatchLevel.NOT_APPLICABLE,
            AlleleMatchLevel.NOT_APPLICABLE
        )

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_levels)

    def test_equal_null_suffix_returns_not_applicable(self):
        """
        Case: two alleles, with risk suffix 'N'
        Allele-1: A*24:09N
        Allele-2: A*24:23N
        Expected: NOT_APPLICABLE
        """
        allele1 = HLA("A*24:09N")
        allele2 = HLA("A*24:09N")
        self.assertEqual(
            allele_match(allele1, allele2),
            AlleleMatchLevel.NOT_APPLICABLE
        )

    def test_identical_g_group_codes_caps_at_ard_match(self):
        """
        Case: both alleles have 'G' group-code
        Allele-1: DQB1*06:02:01G
        Allele-2: DQB1*06:02:02G
        Expected: ARD_MATCH (3)
        """
        a1 = HLA("DQB1*06:02:01G")
        a2 = HLA("DQB1*06:02:02G")
        self.assertEqual(allele_match(a1, a2), AlleleMatchLevel.ARD_MATCH)

    def test_invalid_single_field_p_group_raises(self):
        """
        Case: 'P' group with only one numeric field
        Allele: DQA1*05P
        Expected: MalformedHLAStringError during parsing.
        """
        with self.assertRaises(MalformedHLAStringError):
            HLA("DQA1*05P")

    # TODO: this test is currently not satisfied and will be fixed in the
    # future
    # def test_compressed_allele_without_colon_is_rejected(self):
    #     """
    #     Case: historical 'compressed' notation without ':'
    #     Allele: A*0101
    #     Expected: MalformedHLAStringError.
    #     """
    #     with self.assertRaises(MalformedHLAStringError):
    #         HLA("A*0101")

    def test_non_breaking_space_is_rejected(self):
        """
        Case: allele string contains a non-breaking space (NBSP, U+00A0)
        Allele: 'A*01:01\u00A0'
        Expected: MalformedHLAStringError.
        """
        bad = "A*01:01\u00A0"
        with self.assertRaises(MalformedHLAStringError):
            HLA(bad)

    def test_allele_pair_match_propagates_exceptions(self):
        """
        Test Case: Exception propagation in allele_pair_match
        """
        # Test 1: MalformedHLAStringError
        patient_allele2 = HLA("A*02:01:01")
        donor_allele1 = HLA("B*07:02:01")
        donor_allele2 = HLA("A*02:01:01")

        with self.assertRaises(MalformedHLAStringError):
            # this should already throw error
            patient = HLAPair(hla1=HLA("a"), hla2=patient_allele2)
            donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)
            allele_pair_match(patient=patient, donor=donor)

        # Test 2: InvalidAlleleError
        with self.assertRaises(InvalidAlleleError):
            patient = HLAPair(hla1=HLA("A*07:01"), hla2=patient_allele2)
            donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)
            allele_pair_match(patient=patient, donor=donor)

        # Test 3: TypeError
        patient_allele1 = 1  # Not an HLA object
        with self.assertRaises(TypeError):
            patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
            donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)
            allele_pair_match(patient=patient, donor=donor)

        # Test 4: InvalidLocusComparisonError
        patient_allele1 = HLA("A*01:01:01")
        patient_allele2 = HLA("A*02:01:01")
        with self.assertRaises(InvalidLocusComparisonError):
            patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
            donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)
            allele_pair_match(patient=patient, donor=donor)


# need a correct MatchResult object to test the _loci_level_match
dummy_MatchResult = MatchResult(
                    patient=HLAPair(hla1=HLA('A*01:01'), hla2=HLA('A*01:01')),
                    donor=HLAPair(hla1=HLA('A*01:01'), hla2=HLA('A*01:01')),
                    pairing_score=0,
                    allele_match_levels=[
                        AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH
                    ]
                )


class TestLociLevelMatch_basic_resolution(unittest.TestCase):

    def test_ARD_MATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_ARD_MATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_ARD_MATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_ARD_MATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "PARTIAL_ARD_MISMATCH"
        )

    def test_ARD_MATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "PARTIAL_ARD_MISMATCH"
        )

    def test_ARD_MATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "PARTIAL_ARD_MISMATCH"
        )

    def test_SYNONYMOUS_VARIANT_MATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_NON_CODING_VARIANT_MATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_NON_CODING_VARIANT_MATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_NON_CODING_VARIANT_MATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_NON_CODING_VARIANT_MATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_NON_CODING_VARIANT_MATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_NON_CODING_VARIANT_MATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ALLELE_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ALLELE_MISMATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ALLELE_MISMATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ALLELE_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ALLELE_MISMATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ALLELE_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_TypeError(self):
        level1 = 1
        level2 = 2
        with self.assertRaises(TypeError):
            dummy_MatchResult._calculate_loci_match_basic_resolution(
                level1, level2
            )

    def test_TypeError_message(self):
        level1 = 1
        level2 = 2
        with self.assertRaises(TypeError) as context:
            dummy_MatchResult._calculate_loci_match_basic_resolution(
                level1, level2
            )
        self.assertEqual(
            str(context.exception),
            f"match_level_1 and match_level_2 must be instances of "
            f"{AlleleMatchLevel}, not {type(level1)} and "
            f"{type(level2)}."
        )


class TestLociLevelMatch_high_resolution(unittest.TestCase):

    def test_ARD_MATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_LOCUS_MISMATCH")

    def test_ARD_MATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_GROUP_MISMATCH")

    def test_ARD_MATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_MISMATCH")

    def test_ARD_MATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_ARD_MATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_ARD_MATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_LOCUS_MISMATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_GROUP_MISMATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_MISMATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_SYNONYMOUS_VARIANT_MATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_NON_CODING_VARIANT_MATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_LOCUS_MISMATCH")

    def test_NON_CODING_VARIANT_MATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_GROUP_MISMATCH")

    def test_NON_CODING_VARIANT_MATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_MISMATCH")

    def test_NON_CODING_VARIANT_MATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_NON_CODING_VARIANT_MATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_NON_CODING_VARIANT_MATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MATCH")

    def test_ALLELE_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "LOCUS_MISMATCH_AND_ALLELE_MISMATCH")

    def test_ALLELE_MISMATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ALLELE_GROUP_MISMATCH_AND_ALLELE_MISMATCH")

    def test_ALLELE_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "DOUBLE_ALLELE_MISMATCH")

    def test_ALLELE_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_MISMATCH")

    def test_ALLELE_MISMATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_MISMATCH")

    def test_ALLELE_MISMATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "LOCUS_MISMATCH_AND_ALLELE_GROUP_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "DOUBLE_ALLELE_GROUP_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ALLELE_GROUP_MISMATCH_AND_ALLELE_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_GROUP_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_GROUP_MISMATCH")

    def test_ALLELE_GROUP_MISMATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ALLELE_GROUP_MISMATCH")

    def test_LOCUS_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.LOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "DOUBLE_LOCUS_MISMATCH")

    def test_LOCUS_MISMATCH_and_ALLELE_GROUP_MISMATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "LOCUS_MISMATCH_AND_ALLELE_GROUP_MISMATCH")

    def test_LOCUS_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "LOCUS_MISMATCH_AND_ALLELE_MISMATCH")

    def test_LOCUS_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_LOCUS_MISMATCH")

    def test_LOCUS_MISMATCH_and_SYNONYMOUS_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_LOCUS_MISMATCH")

    def test_LOCUS_MISMATCH_and_NON_CODING_VARIANT_MATCH(self):
        level1 = AlleleMatchLevel.LOCUS_MISMATCH
        level2 = AlleleMatchLevel.NON_CODING_VARIANT_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_LOCUS_MISMATCH")

    def test_TypeError(self):
        level1 = 1
        level2 = 2
        with self.assertRaises(TypeError):
            dummy_MatchResult._calculate_loci_match_basic_resolution(
                level1, level2
            )

    def test_TypeError_message(self):
        level1 = 1
        level2 = 2
        with self.assertRaises(TypeError) as context:
            dummy_MatchResult._calculate_loci_match_basic_resolution(
                level1, level2
            )
        self.assertEqual(
            str(context.exception),
            f"match_level_1 and match_level_2 must be instances of "
            f"{AlleleMatchLevel}, not {type(level1)} and "
            f"{type(level2)}."
        )


class TestMultiLocusMatch(unittest.TestCase):

    def test_multi_locus_match(self):
        """
        Test match results for multiple pairs of alleles stored in an
        individual object
        """
        patient_allele1 = HLA("A*01:01:01")
        patient_allele2 = HLA("A*01:01:01")

        patient_hla_pair1 = HLAPair(hla1=patient_allele1, hla2=patient_allele2)

        patient_allele1 = HLA("B*07:02:01")
        patient_allele2 = HLA("B*07:02:01")

        patient_hla_pair2 = HLAPair(hla1=patient_allele1, hla2=patient_allele2)

        patient = Individual(hla_data=[patient_hla_pair1, patient_hla_pair2])

        donor_allele1 = HLA("A*01:01:01")
        donor_allele2 = HLA("A*01:01:01")

        donor_hla_pair1 = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        donor_allele1 = HLA("B*07:02:01")
        donor_allele2 = HLA("B*07:02:01")

        donor_hla_pair2 = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        donor = Individual(hla_data=[donor_hla_pair1, donor_hla_pair2])

        result = multi_locus_match(patient, donor)

        self.assertEqual(
            result[0].allele_match_levels,
            (
                AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH,
                AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
            )
        )
        self.assertEqual(
            result[1].allele_match_levels,
            (
                AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH,
                AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
            )
        )

    def test_multi_locus_match_missing_donor_locus(self):
        patient_allele1 = HLA("A*01:01:01")
        patient_allele2 = HLA("A*01:01:01")

        patient_hla_pair1 = HLAPair(hla1=patient_allele1, hla2=patient_allele2)

        patient_allele1 = HLA("B*07:02:01")
        patient_allele2 = HLA("B*07:02:01")

        patient_hla_pair2 = HLAPair(hla1=patient_allele1, hla2=patient_allele2)

        patient = Individual(hla_data=[patient_hla_pair1, patient_hla_pair2])

        donor_allele1 = HLA("A*01:01:01")
        donor_allele2 = HLA("A*01:01:01")

        donor_hla_pair1 = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        # different locus recorded for donor
        donor_allele1 = HLA("C*07:02:01")
        donor_allele2 = HLA("C*07:02:01")

        donor_hla_pair2 = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        donor = Individual(hla_data=[donor_hla_pair1, donor_hla_pair2])

        result = multi_locus_match(patient, donor)

        # check if correct warning were raised during execution
        with self.assertLogs("py_hla_match.matching", level="WARNING") as cm:
            result = multi_locus_match(patient, donor)
            self.assertEqual(len(cm.output), 2)
            self.assertIn(
                "matching will be reported as NOT_APPLICABLE",
                cm.output[0],
            )

        self.assertEqual(len(result), 2)
        self.assertEqual(
            result[0].allele_match_levels,
            (
                AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH,
                AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
            )
        )

    def test_donor_extra_locus_is_ignored(self):
        # patient: locus A only
        patient_pair = HLAPair(HLA("A*01:01"), HLA("A*01:01"))
        patient = Individual([patient_pair])

        # donor: locus A + extra locus B
        donor_pair_a = HLAPair(HLA("A*01:01"), HLA("A*01:01"))
        donor_pair_b = HLAPair(HLA("B*07:02"), HLA("B*07:02"))
        donor = Individual([donor_pair_a, donor_pair_b])

        result = multi_locus_match(patient, donor)
        self.assertEqual(len(result), 1)            # only locus A evaluated
        self.assertEqual(result[0].patient.locus, "A")

    def test_insufficient_resolution_logs_warning(self):
        patient_pair = HLAPair(HLA("C*NA"), HLA("C*NA"))
        donor_pair = HLAPair(HLA("C*01:02"), HLA("C*01:02"))
        patient = Individual([patient_pair])
        donor = Individual([donor_pair])

        with self.assertLogs("py_hla_match.matching", level="WARNING") as cm:
            result = multi_locus_match(patient, donor)

        # at least one warning about insufficient resolution
        self.assertTrue(
            any("insufficient" in msg.lower() for msg in cm.output)
        )
        self.assertEqual(
            result[0].allele_match_levels,
            (
                AlleleMatchLevel.NOT_APPLICABLE,
                AlleleMatchLevel.NOT_APPLICABLE
            ),
        )


if __name__ == "__main__":
    unittest.main()
