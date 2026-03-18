import unittest

from py_hla_match.hla import HLA
from py_hla_match.matching import (
    allele_match,
    _get_correct_allele_pairing,
    allele_pair_match,
    MatchResult,
    multi_locus_match
)
from py_hla_match.exceptions import (
    InvalidLocusComparisonError, MalformedHLAStringError
)
from pyard.exceptions import InvalidAlleleError
from py_hla_match.models import HLAPair, Individual
from py_hla_match.policy import (
    AlleleMatchLevel,
    ARDMatchLevel,
    ARDMatchLevelCertainty,
    MolecularMatchLevel,
    MolecularMatchLevelCertainty,
    ExpressionSuffixPolicy,
    ExpressionSuffixMatchLevel
)
from py_hla_match.config import (
    HLAMatchConfig,
    set_config
)


class TestAlleleMatch(unittest.TestCase):

    def tearDown(self) -> None:
        set_config(HLAMatchConfig())

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
        Expected Match Level: DRB345_SUBLOCUS_MISMATCH
        """
        allele1 = HLA("DRB3*02:02:01")
        allele2 = HLA("DRB4*01:03:01")
        expected_match_level = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_drb35_locus_mismatch(self):
        """
        Test Case: Locus mismatch
        Allele 1: DRB3*01
        Allele 2: DRB5*01
        Expected Match Level: DRB345_SUBLOCUS_MISMATCH
        """
        allele1 = HLA("DRB3*01")
        allele2 = HLA("DRB5*01")
        expected_match_level = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_drb3x_locus_mismatch(self):
        """
        Test Case: Locus mismatch
        Allele 1: DRB3*01
        Allele 2: DRBX*NE
        Expected Match Level: DRB345_SUBLOCUS_MISMATCH
        """
        allele1 = HLA("DRB3*01")
        allele2 = HLA("DRBX*NE")
        expected_match_level = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
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
        Expected Match Level: ANTIGEN_MISMATCH
        """
        allele1 = HLA("DPB1*02:01:02")
        allele2 = HLA("DPB1*04:02:01")
        expected_match_level = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_allele_mismatch_with_P_group(self):
        """
        Test Case: Allele mismatch with 'P' suffix
        Allele 1: DPB1*04:01P
        Allele 2: DPB1*04:02P
        Expected Match Level: ALLELE_MISMATCH
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
        Expected Match Level: ARD_MATCH
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
        Expected Match Level: ARD_MATCH
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
        Expected Match Level: NOT_ASSESSABLE
        """
        allele1 = HLA("A*01:436Q")
        allele2 = HLA("A*01:01:70")
        expected_match_level = AlleleMatchLevel.NOT_ASSESSABLE
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_valid_allele_mismatch_with_L_suffix(self):
        """
        Test Case: ARD match with L suffix
        Allele 1: B*38:68L
        Allele 2: B*38:01P
        Expected Match Level: ALLELE_MISMATCH
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
        Expected Match Level: ALLELE_MISMATCH
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
        Expected Match Level: ARD_MATCH
        """
        allele1 = HLA("DPA1*01:03:01")
        allele2 = HLA("DPA1*01:03:01")
        expected_match_level = AlleleMatchLevel.ARD_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_cap_at_ard_match(self):
        """
        Test Case: Non encoding variant match
        Allele 1: A*01:01:01:46
        Allele 2: A*01:01:01:46
        Expected Match Level: ARD_MATCH
        """
        allele1 = HLA("A*01:01:01:46")
        allele2 = HLA("A*01:01:01:46")
        expected_match_level = AlleleMatchLevel.ARD_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_not_applicable_due_to_nan(self):
        """
        Test Case: missing allele data *NE
        Allele 1: A*NE
        Allele 2: A*01:01
        Expected Match Level: NOT_ASSESSABLE
        """
        allele1 = HLA("A*NE")
        allele2 = HLA("A*01:01")
        expected_match_level = AlleleMatchLevel.NOT_ASSESSABLE
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_one_field_vs_two_field_same_group_not_applicable(self):
        """
        Test Case: 1-field vs 2-field with same group
        Allele 1: B*07
        Allele 2: B*07:05
        Expected Match Level: NOT_ASSESSABLE
        """
        allele1 = HLA("B*07")
        allele2 = HLA("B*07:05")
        expected_match_level = AlleleMatchLevel.NOT_ASSESSABLE
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_allele_group_mismatch_with_one_field(self):
        """
        Test Case: 1-field vs 1-field with different groups
        Allele 1: C*01
        Allele 2: C*02
        Expected Match Level: ANTIGEN_MISMATCH
        """
        allele1 = HLA("C*01")
        allele2 = HLA("C*02")
        expected_match_level = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_not_applicable_with_suffix_Q(self):
        """
        Test Case: Q suffix default is NOT_ASSESSABLE
        Allele 1: A*24:473Q
        Allele 2: A*02:99
        Expected Match Level: NOT_ASSESSABLE
        """
        allele1 = HLA("A*24:473Q")
        allele2 = HLA("A*24:02P")
        expected_match_level = AlleleMatchLevel.NOT_ASSESSABLE
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

    def test_synonymous_vs_full_variant_cap_at_ard(self):
        """
        Test Case: 3-field vs 4-field with equal synonymous variant
        Allele 1: DQA1*01:01:02
        Allele 2: DQA1*01:01:02:07
        Expected Match Level: ARD_MATCH
        """
        allele1 = HLA("DQA1*01:01:02")
        allele2 = HLA("DQA1*01:01:02:05")
        expected_match_level = AlleleMatchLevel.ARD_MATCH
        result = allele_match(allele1, allele2)
        self.assertEqual(result, expected_match_level)

    def test_equal_risk_suffix_policy_can_be_set_to_mismatch(self):
        # Configure equal risk (e.g., N vs N) to ALLELE_MISMATCH instead of
        # default NOT_ASSESSABLE
        pol = ExpressionSuffixPolicy(
            equal_risk=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            risk_vs_none=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            risk_vs_different_risk=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            q_present=ExpressionSuffixMatchLevel.NOT_ASSESSABLE,
        )
        set_config(HLAMatchConfig(expression_suffix_policy=pol))
        a1 = HLA("A*24:09N")
        a2 = HLA("A*24:09N")
        self.assertEqual(
            allele_match(a1, a2), AlleleMatchLevel.ALLELE_MISMATCH
        )

    def test_q_present_policy_can_be_set_to_mismatch(self):
        # Configure Q present to be treated as mismatch rather than default
        # NOT_ASSESSABLE
        pol = ExpressionSuffixPolicy(
            equal_risk=ExpressionSuffixMatchLevel.NOT_ASSESSABLE,
            risk_vs_none=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            risk_vs_different_risk=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            q_present=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
        )
        set_config(HLAMatchConfig(expression_suffix_policy=pol))
        a1 = HLA("A*01:436Q")
        a2 = HLA("A*01:01:70")
        self.assertEqual(
            allele_match(a1, a2), AlleleMatchLevel.ALLELE_MISMATCH
        )

    def test_risk_vs_none_policy_can_be_set_to_not_applicable(self):
        # Configure risk vs none (e.g., N vs none) to NOT_ASSESSABLE after
        # ARD-equivalence
        pol = ExpressionSuffixPolicy(
            equal_risk=ExpressionSuffixMatchLevel.NOT_ASSESSABLE,
            risk_vs_none=ExpressionSuffixMatchLevel.NOT_ASSESSABLE,
            risk_vs_different_risk=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            q_present=ExpressionSuffixMatchLevel.NOT_ASSESSABLE,
        )
        set_config(HLAMatchConfig(expression_suffix_policy=pol))
        # ARD-equivalent: A*24:09N reduces to A*24:02; compared to A*24:02
        a1 = HLA("A*24:09N")
        a2 = HLA("A*24:02")
        self.assertEqual(allele_match(a1, a2), AlleleMatchLevel.NOT_ASSESSABLE)

    def test_suffix_policy_does_not_override_group_mismatch_severity(self):
        # Even with aggressive suffix policy, group mismatch must dominate
        pol = ExpressionSuffixPolicy(
            equal_risk=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            risk_vs_none=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            risk_vs_different_risk=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
            q_present=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
        )
        set_config(HLAMatchConfig(expression_suffix_policy=pol))
        a1 = HLA("A*23:41")   # group 23
        a2 = HLA("A*24:09N")   # group 24
        self.assertEqual(
            allele_match(a1, a2), AlleleMatchLevel.ANTIGEN_MISMATCH
        )


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

        result = _get_correct_allele_pairing(patient, donor)

        expected_pairing = (
            allele_match(p1, d1),
            allele_match(p2, d2),
        )
        self.assertEqual(result.allele_match_levels, expected_pairing)
        self.assertEqual(result.score, sum(expected_pairing))

    def test_all_not_applicable_pairing_score(self):
        """
        Test Case: NOT_ASSESSABLE twice
        Patient Alleles: B*NA, B*NA
        Donor Alleles:   B*NA, B*NA
        Expected Best Score: AlleleMatchLevel.NOT_ASSESSABLE * 2
        Expected Pairing Levels:
            (AlleleMatchLevel.NOT_ASSESSABLE, AlleleMatchLevel.NOT_ASSESSABLE)
        """
        p1 = p2 = HLA("B*NA")
        d1 = d2 = HLA("B*NA")
        patient = HLAPair(p1, p2)
        donor = HLAPair(d1, d2)

        result = _get_correct_allele_pairing(patient, donor)

        self.assertEqual(result.score, AlleleMatchLevel.NOT_ASSESSABLE * 2)
        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.NOT_ASSESSABLE, AlleleMatchLevel.NOT_ASSESSABLE),
        )

    def test_pairing_prefers_lower_negative_penalty(self):
        """
        Test Case: Two possible pairings – one sums to 0, the other to +1.
        """
        p = HLAPair(HLA("B*07"),     HLA("B*07:02"))
        d = HLAPair(HLA("B*07:02"),  HLA("B*07"))

        result = _get_correct_allele_pairing(p, d)
        self.assertEqual(result.score, 1)
        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.NOT_ASSESSABLE, AlleleMatchLevel.ARD_MATCH)
        )

    def test_correct_pairing_with_ambiguous_alleles(self):
        """
        Test Case: pairing prefers NOT_ASSESSABLE over any mismatch
        """
        p = HLAPair(HLA("B*07"),     HLA("B*01"))
        d = HLAPair(HLA("B*07"),  HLA("B*01"))

        result = _get_correct_allele_pairing(p, d)

        self.assertEqual(result.score, 0)
        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.NOT_ASSESSABLE, AlleleMatchLevel.NOT_ASSESSABLE)
        )


class TestMatchResultARDRefinementFieldsDefaults(unittest.TestCase):
    """Tests for ARD/molecular fields in MatchResult construction."""

    def test_defaults_when_molecular_fields_not_provided(self):
        """
        If MatchResult is constructed with only allele_match_levels:

        - ard_match_levels / ard_match_certainties
        - molecular_match_levels / molecular_match_certainties

        must all default to NOT_APPLICABLE. They are only meaningful when
        computed via allele_pair_match (i.e. AlleleMatchLevel == ARD_MATCH).
        """
        patient = HLAPair(HLA("A*01:01"), HLA("A*01:01"))
        donor = HLAPair(HLA("A*01:01"), HLA("A*01:01"))

        result = MatchResult(
            patient=patient,
            donor=donor,
            pairing_score=AlleleMatchLevel.ARD_MATCH * 2,
            allele_match_levels=(
                AlleleMatchLevel.ARD_MATCH,
                AlleleMatchLevel.ARD_MATCH,
            ),
        )

        self.assertEqual(
            result.ard_match_levels,
            (
                ARDMatchLevel.NOT_APPLICABLE,
                ARDMatchLevel.NOT_APPLICABLE
            ),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (
                ARDMatchLevelCertainty.NOT_APPLICABLE,
                ARDMatchLevelCertainty.NOT_APPLICABLE
            ),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (
                MolecularMatchLevel.NOT_APPLICABLE,
                MolecularMatchLevel.NOT_APPLICABLE
            ),
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (
                MolecularMatchLevelCertainty.NOT_APPLICABLE,
                MolecularMatchLevelCertainty.NOT_APPLICABLE
            ),
        )

    def test_explicit_molecular_fields_are_preserved(self):
        """
        If ARD/molecular fields are passed explicitly to MatchResult,
        they must be stored unchanged and NOT recomputed.
        """
        patient = HLAPair(HLA("A*01:01:01:04"), HLA("A*01:01:01:01"))
        donor = HLAPair(HLA("A*01:01:01:01"), HLA("A*01:01:01:03"))

        result = MatchResult(
            patient=patient,
            donor=donor,
            pairing_score=AlleleMatchLevel.ARD_MATCH * 2,
            allele_match_levels=(
                AlleleMatchLevel.ARD_MATCH,
                AlleleMatchLevel.ARD_MATCH,
            ),
            ard_match_levels=(
                ARDMatchLevel.G_GROUP_MATCH,
                ARDMatchLevel.G_GROUP_MATCH,
            ),
            ard_match_level_certainty=(
                ARDMatchLevelCertainty.CERTAIN,
                ARDMatchLevelCertainty.CERTAIN,
            ),
            molecular_match_levels=(
                MolecularMatchLevel.EXACT_ALLELE_MATCH,
                MolecularMatchLevel.CODING_SEQUENCE_MATCH,
            ),
            molecular_match_level_certainty=(
                MolecularMatchLevelCertainty.CERTAIN,
                MolecularMatchLevelCertainty.CERTAIN,
            ),
        )

        self.assertEqual(
            result.ard_match_levels,
            (
                ARDMatchLevel.G_GROUP_MATCH,
                ARDMatchLevel.G_GROUP_MATCH
            ),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (
                ARDMatchLevelCertainty.CERTAIN,
                ARDMatchLevelCertainty.CERTAIN
            ),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (
                MolecularMatchLevel.EXACT_ALLELE_MATCH,
                MolecularMatchLevel.CODING_SEQUENCE_MATCH
            ),
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (
                MolecularMatchLevelCertainty.CERTAIN,
                MolecularMatchLevelCertainty.CERTAIN
            ),
        )


class TestAllelePairMatch(unittest.TestCase):
    def test_valid_match_without_swapping(self):
        """
        Test Case: Double ARD_MATCH (4)
        Patient: DRB1*15:01:01, DRB1*15:01:01
        Donor: DRB1*15:01:01, DRB1*15:01:01
        Expected Score: ARD_MATCH + ARD_MATCH
        Expected Allele Match Levels: Double ARD_MATCH
        """
        patient_allele1 = HLA("DRB1*15:01:01")
        patient_allele2 = HLA("DRB1*15:01:01")
        donor_allele1 = HLA("DRB1*15:01:01")
        donor_allele2 = HLA("DRB1*15:01:01")

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

    def test_valid_match_with_suffixes_G_and_P(self):
        """
        Test Case: Match with suffixes 'G' and 'P'
        Patient: DRB1*01:01:01G, DRB1*07:01:01G
        Donor: DRB1*01:01P, DRB1*07:01:01
        Expected Score: ARD_MATCH + ARD_MATCH
        Expected Allele Match Levels: Double ARD_MATCH
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
        Expected Score: ARD_MATCH + ARD_MATCH
        Expected Allele Match Levels:
            - B*35:02:01 with B*35:02:01: ARD_MATCH
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
            AlleleMatchLevel.ARD_MATCH,
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
        Test Case: Double ALLELE_MISMATCH
        Patient Alleles: DPB1*04:01:01, DPB1*04:01:01
        Donor Alleles: DPB1*04:02:01, DPB1*04:02:01
        Expected Score: ALLELE_MISMATCH * 2
        Expected Allele Match Levels: Double ALLELE_MISMATCH
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
        Test Case: Swapping ANTIGEN_MISMATCH, ALLELE_MISMATCH
        Patient Alleles: DPB1*01:01:01, DPB1*04:02:01
        Donor Alleles: DPB1*04:01:01, DPB1*02:01:02
        Expected Score: ALLELE_MISMATCH + ANTIGEN_MISMATCH
        Expected Allele Match Levels: [ALLELE_MISMATCH, ANTIGEN_MISMATCH]
        """
        patient_allele1 = HLA("DPB1*01:01:01")
        patient_allele2 = HLA("DPB1*04:02:01")
        donor_allele1 = HLA("DPB1*04:01:01")
        donor_allele2 = HLA("DPB1*02:01:02")

        patient = HLAPair(hla1=patient_allele1, hla2=patient_allele2)
        donor = HLAPair(hla1=donor_allele1, hla2=donor_allele2)

        result = allele_pair_match(patient, donor)

        expected_match_levels = (
            AlleleMatchLevel.ANTIGEN_MISMATCH,
            AlleleMatchLevel.ALLELE_MISMATCH
        )
        expected_score = sum(expected_match_levels)

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_match_levels)

    def test_valid_ard_match_and_allele_mismatch_swapping(self):
        """
        Test Case: ARD_MATCH and ALLELE_MISMATCH with swapping
        Patient Alleles: DPB1*03:01P, DPB1*04:01P
        Donor Alleles: 5 DPB1*04:02P, DPB1*03:01P
        Expected Score: ARD_MATCH + ALLELE_MISMATCH
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
        Test Case: Both comparisons NOT_ASSESSABLE
        Patient Alleles: E*NE, E*NE
        Donor Alleles:   E*NE, E*NE
        Expected Pairing Score: AlleleMatchLevel.NOT_ASSESSABLE * 2
        Expected Match Levels: (NOT_ASSESSABLE, NOT_ASSESSABLE)
        """
        na1 = HLA("E*NE")
        na2 = HLA("E*NE")
        na3 = HLA("E*NE")
        na4 = HLA("E*NE")
        patient = HLAPair(na1, na2)
        donor = HLAPair(na3, na4)

        result = allele_pair_match(patient, donor)

        expected_score = AlleleMatchLevel.NOT_ASSESSABLE * 2
        expected_levels = (
            AlleleMatchLevel.NOT_ASSESSABLE,
            AlleleMatchLevel.NOT_ASSESSABLE
        )

        self.assertEqual(result.pairing_score, expected_score)
        self.assertEqual(result.allele_match_levels, expected_levels)

    def test_equal_null_suffix_returns_not_applicable(self):
        """
        Case: two alleles, with risk suffix 'N'
        Allele-1: A*24:09N
        Allele-2: A*24:23N
        Expected: NOT_ASSESSABLE
        """
        allele1 = HLA("A*24:09N")
        allele2 = HLA("A*24:09N")
        self.assertEqual(
            allele_match(allele1, allele2),
            AlleleMatchLevel.NOT_ASSESSABLE
        )

    def test_identical_g_group_codes_caps_at_ard_match(self):
        """
        Case: both alleles have 'G' group-code
        Allele-1: DQB1*06:02:01G
        Allele-2: DQB1*06:02:01G
        Expected: ARD_MATCH
        """
        a1 = HLA("DQB1*06:02:01G")
        a2 = HLA("DQB1*06:02:01G")
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


class TestAllelePairMatchARDRefinement(unittest.TestCase):
    """
    Tests for ARD refinement including both MolecularMatchLevel and
    ARDMatchLevel computation.

    Tests validate refinement for cases AlleleMatchLevel is ARD_MATCH and
    ensure that non-ARD_MATCH cases map to NOT_APPLICABLE.
    """

    def tearDown(self) -> None:
        set_config(HLAMatchConfig())

    # NOT_APPLICABLE cases (AlleleMatchLevel != ARD_MATCH)
    def test_antigen_mismatch_molecular_not_applicable(self):
        """
        ANTIGEN_MISMATCH:
        - ARDMatchLevel / ARDMatchLevelCertainty -> NOT_APPLICABLE
        - MolecularMatchLevel / MolecularMatchLevelCertainty -> NOT_APPLICABLE
        """
        patient = HLAPair(HLA("C*01:02"), HLA("C*01:12"))
        donor = HLAPair(HLA("C*03:02:01"), HLA("C*02:02"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ANTIGEN_MISMATCH,
             AlleleMatchLevel.ANTIGEN_MISMATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.NOT_APPLICABLE,
             ARDMatchLevel.NOT_APPLICABLE),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.NOT_APPLICABLE,
             ARDMatchLevelCertainty.NOT_APPLICABLE),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.NOT_APPLICABLE,
             MolecularMatchLevel.NOT_APPLICABLE)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.NOT_APPLICABLE,
             MolecularMatchLevelCertainty.NOT_APPLICABLE)
        )

    def test_allele_mismatch_molecular_not_applicable(self):
        """
        ALLELE_MISMATCH:
        - ARD refinement is not applied
        - Molecular refinement is not applied
        => all ARD/molecular fields NOT_APPLICABLE
        """
        patient = HLAPair(HLA("DPB1*04:01:01:01"), HLA("DPB1*02:02:01"))
        donor = HLAPair(HLA("DPB1*02:01:02:01"), HLA("DPB1*04:02:01"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ALLELE_MISMATCH,
             AlleleMatchLevel.ALLELE_MISMATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.NOT_APPLICABLE,
             ARDMatchLevel.NOT_APPLICABLE),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.NOT_APPLICABLE,
             ARDMatchLevelCertainty.NOT_APPLICABLE),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.NOT_APPLICABLE,
             MolecularMatchLevel.NOT_APPLICABLE)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.NOT_APPLICABLE,
             MolecularMatchLevelCertainty.NOT_APPLICABLE)
        )

    def test_not_assessable_molecular_not_applicable(self):
        """
        NOT_ASSESSABLE (insufficient typing resolution):
        - ARD refinement is not applied
        - Molecular refinement is not applied
        => all ARD/molecular fields NOT_APPLICABLE
        """
        patient = HLAPair(HLA("B*07"), HLA("B*08"))
        donor = HLAPair(HLA("B*08:01"), HLA("B*07:02"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.NOT_ASSESSABLE,
             AlleleMatchLevel.NOT_ASSESSABLE)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.NOT_APPLICABLE,
             ARDMatchLevel.NOT_APPLICABLE),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.NOT_APPLICABLE,
             ARDMatchLevelCertainty.NOT_APPLICABLE),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.NOT_APPLICABLE,
             MolecularMatchLevel.NOT_APPLICABLE)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.NOT_APPLICABLE,
             MolecularMatchLevelCertainty.NOT_APPLICABLE)
        )

    # EXACT_ALLELE_MATCH cases
    def test_exact_allele_match_four_field_identical(self):
        """
        Four-field identical alleles (A*01:01:01:01 vs A*01:01:01:01):
        - ARD refinement: G_GROUP_MATCH, CERTAIN
        - Molecular: EXACT_ALLELE_MATCH, CERTAIN

        """
        patient = HLAPair(HLA("A*01:01:01:01"), HLA("A*02:01:01:01"))
        donor = HLAPair(HLA("A*02:01:01:01"), HLA("A*01:01:01:01"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.G_GROUP_MATCH,
             ARDMatchLevel.G_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.CERTAIN,
             ARDMatchLevelCertainty.CERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.EXACT_ALLELE_MATCH,
             MolecularMatchLevel.EXACT_ALLELE_MATCH)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.CERTAIN,
             MolecularMatchLevelCertainty.CERTAIN)
        )

    # CODING_SEQUENCE_MATCH cases
    def test_coding_sequence_match_fourth_field_differs(self):
        """
        Fourth field differs, same 1–3 fields:
        - ARD refinement: G_GROUP_MATCH, CERTAIN
        - Molecular: CODING_SEQUENCE_MATCH, CERTAIN
        """
        patient = HLAPair(HLA("A*01:01:01:01"), HLA("A*01:01:01:04"))
        donor = HLAPair(HLA("A*01:01:01:03"), HLA("A*01:01:01:05"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.G_GROUP_MATCH,
             ARDMatchLevel.G_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.CERTAIN,
             ARDMatchLevelCertainty.CERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.CODING_SEQUENCE_MATCH,
             MolecularMatchLevel.CODING_SEQUENCE_MATCH)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.CERTAIN,
             MolecularMatchLevelCertainty.CERTAIN)
        )

    def test_coding_sequence_match_three_field_vs_four_field_uncertain(self):
        """
        Three-field vs four-field, same 1–3 fields:
        - ARD refinement: G_GROUP_MATCH, CERTAIN
        - Molecular: CODING_SEQUENCE_MATCH, UNCERTAIN
        """
        patient = HLAPair(HLA("A*01:01:01"), HLA("A*01:02:01"))
        donor = HLAPair(HLA("A*01:02:01:01"), HLA("A*01:01:01:03"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.G_GROUP_MATCH,
             ARDMatchLevel.G_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.CERTAIN,
             ARDMatchLevelCertainty.CERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.CODING_SEQUENCE_MATCH,
             MolecularMatchLevel.CODING_SEQUENCE_MATCH)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.UNCERTAIN,
             MolecularMatchLevelCertainty.UNCERTAIN)
        )

    # FULL_PROTEIN_MATCH cases
    def test_full_protein_match_third_field_differs(self):
        """
        Third field differs, same 1–2 fields:
        - ARD refinement: P_GROUP_MATCH, UNCERTAIN
        - Molecular: FULL_PROTEIN_MATCH, CERTAIN
        """
        patient = HLAPair(HLA("A*01:01:02"), HLA("A*01:01:01"))
        donor = HLAPair(HLA("A*01:01:03"), HLA("A*01:01:04"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.P_GROUP_MATCH,
             ARDMatchLevel.P_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.UNCERTAIN,
             ARDMatchLevelCertainty.UNCERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.FULL_PROTEIN_MATCH,
             MolecularMatchLevel.FULL_PROTEIN_MATCH)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.CERTAIN,
             MolecularMatchLevelCertainty.CERTAIN)
        )

    def test_full_protein_match_two_field_vs_three_field_uncertain(self):
        """
        Two-field vs three-field with identical 1–2 fields:
        - ARD refinement: P_GROUP_MATCH, UNCERTAIN
        - Molecular: FULL_PROTEIN_MATCH, UNCERTAIN
        """
        patient = HLAPair(HLA("A*01:02"), HLA("A*01:01"))
        donor = HLAPair(HLA("A*01:01:01"), HLA("A*01:02:01"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.P_GROUP_MATCH,
             ARDMatchLevel.P_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.UNCERTAIN,
             ARDMatchLevelCertainty.UNCERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.FULL_PROTEIN_MATCH,
             MolecularMatchLevel.FULL_PROTEIN_MATCH)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.UNCERTAIN,
             MolecularMatchLevelCertainty.UNCERTAIN)
        )

    def test_full_protein_match_two_field_identical_uncertain(self):
        """
        Two-field identical alleles:
        - ARD refinement: P_GROUP_MATCH, UNCERTAIN
        - Molecular: FULL_PROTEIN_MATCH, UNCERTAIN

        """
        patient = HLAPair(HLA("A*01:01"), HLA("A*01:02"))
        donor = HLAPair(HLA("A*01:02"), HLA("A*01:01"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.P_GROUP_MATCH,
             ARDMatchLevel.P_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.UNCERTAIN,
             ARDMatchLevelCertainty.UNCERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.FULL_PROTEIN_MATCH,
             MolecularMatchLevel.FULL_PROTEIN_MATCH)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.UNCERTAIN,
             MolecularMatchLevelCertainty.UNCERTAIN)
        )

    # NOT_ASSESSABLE cases (P-group)
    def test_p_group_molecular_not_assessable(self):
        """
        P-group alleles present (A*01:02P, A*01:01P):
        - ARD refinement: P_GROUP_MATCH, UNCERTAIN
        - Molecular: NOT_ASSESSABLE, UNCERTAIN
        """
        patient = HLAPair(HLA("A*01:02P"), HLA("A*01:01P"))
        donor = HLAPair(HLA("A*01:01:01"), HLA("A*01:412"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.P_GROUP_MATCH,
             ARDMatchLevel.P_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.UNCERTAIN,
             ARDMatchLevelCertainty.UNCERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.NOT_ASSESSABLE,
             MolecularMatchLevel.NOT_ASSESSABLE)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.UNCERTAIN,
             MolecularMatchLevelCertainty.UNCERTAIN)
        )

    # G-group cases
    def test_g_group_vs_specific_full_protein_match(self):
        """
        G-group vs specific allele within same G-group:
        - ARD refinement: G_GROUP_MATCH, CERTAIN
        - Molecular: NOT_ASSESSABLE, UNCERTAIN
        """
        patient = HLAPair(HLA("C*07:02:01G"), HLA("C*07:02:10G"))
        donor = HLAPair(HLA("C*07:02:10"), HLA("C*07:02:01"))
        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH)
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.G_GROUP_MATCH,
             ARDMatchLevel.G_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (ARDMatchLevelCertainty.CERTAIN,
             ARDMatchLevelCertainty.CERTAIN),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (MolecularMatchLevel.NOT_ASSESSABLE,
             MolecularMatchLevel.NOT_ASSESSABLE)
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (MolecularMatchLevelCertainty.UNCERTAIN,
             MolecularMatchLevelCertainty.UNCERTAIN)
        )

    def test_a0201_p_group_two_field_ard_match_only(self):
        """
        Two-field alleles A*02:66 and A*02:75 in A*02:01P:
        - ARD refinement: P_GROUP_MATCH, UNCERTAIN
        - Molecular: ARD_MATCH_ONLY, CERTAIN
        """
        patient = HLAPair(HLA("A*02:66"), HLA("A*02:66"))
        donor = HLAPair(HLA("A*02:75"), HLA("A*02:75"))

        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH),
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.P_GROUP_MATCH, ARDMatchLevel.P_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (
                ARDMatchLevelCertainty.UNCERTAIN,
                ARDMatchLevelCertainty.UNCERTAIN,
            ),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (
                MolecularMatchLevel.ARD_MATCH_ONLY,
                MolecularMatchLevel.ARD_MATCH_ONLY,
            ),
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (
                MolecularMatchLevelCertainty.CERTAIN,
                MolecularMatchLevelCertainty.CERTAIN,
            ),
        )

    def test_ard_match_only_second_field_differs_with_three_field(self):
        """
        Two 3-field alleles in the same ARD group but with different 2nd field
        - AlleleMatchLevel: ARD_MATCH
        - ARD refinement: G_GROUP_MATCH, CERTAIN
        - Molecular: ARD_MATCH_ONLY, CERTAIN
        """
        patient = HLAPair(HLA("A*02:01:01"), HLA("A*02:01:01"))
        donor = HLAPair(HLA("A*02:09:01"), HLA("A*02:09:01"))

        result = allele_pair_match(patient, donor)

        self.assertEqual(
            result.allele_match_levels,
            (AlleleMatchLevel.ARD_MATCH, AlleleMatchLevel.ARD_MATCH),
        )
        self.assertEqual(
            result.ard_match_levels,
            (ARDMatchLevel.G_GROUP_MATCH, ARDMatchLevel.G_GROUP_MATCH),
        )
        self.assertEqual(
            result.ard_match_certainties,
            (
                ARDMatchLevelCertainty.CERTAIN,
                ARDMatchLevelCertainty.CERTAIN,
            ),
        )
        self.assertEqual(
            result.molecular_match_levels,
            (
                MolecularMatchLevel.ARD_MATCH_ONLY,
                MolecularMatchLevel.ARD_MATCH_ONLY,
            ),
        )
        self.assertEqual(
            result.molecular_match_certainties,
            (
                MolecularMatchLevelCertainty.CERTAIN,
                MolecularMatchLevelCertainty.CERTAIN,
            ),
        )


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

    def test_ARD_MATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "PARTIAL_ARD_MISMATCH"
        )

    def test_ARD_MATCH_and_ANTIGEN_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "PARTIAL_ARD_MISMATCH"
        )

    def test_ARD_MATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ARD_MATCH
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "PARTIAL_ARD_MISMATCH"
        )

    def test_ALLELE_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
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

    def test_ALLELE_MISMATCH_and_ANTIGEN_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ALLELE_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ANTIGEN_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_ANTIGEN_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ANTIGEN_MISMATCH_and_ANTIGEN_MISMATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_ANTIGEN_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_ANTIGEN_MISMATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_basic_resolution(
            level1, level2
        )
        self.assertEqual(result, "ARD_MISMATCH")

    def test_LOCUS_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
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
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_DRB345_SUBLOCUS_MISMATCH")

    def test_ALLELE_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "DRB345_SUBLOCUS_MISMATCH_AND_ALLELE_MISMATCH"
        )

    def test_ALLELE_MISMATCH_and_ANTIGEN_MISMATCH(self):
        level1 = AlleleMatchLevel.ALLELE_MISMATCH
        level2 = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ANTIGEN_MISMATCH_AND_ALLELE_MISMATCH")

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

    def test_ANTIGEN_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "DRB345_SUBLOCUS_MISMATCH_AND_ANTIGEN_MISMATCH"
        )

    def test_ANTIGEN_MISMATCH_and_ANTIGEN_MISMATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "DOUBLE_ANTIGEN_MISMATCH")

    def test_ANTIGEN_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "ANTIGEN_MISMATCH_AND_ALLELE_MISMATCH")

    def test_ANTIGEN_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.ANTIGEN_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_ANTIGEN_MISMATCH")

    def test_LOCUS_MISMATCH_and_LOCUS_MISMATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "DOUBLE_DRB345_SUBLOCUS_MISMATCH")

    def test_LOCUS_MISMATCH_and_ANTIGEN_MISMATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.ANTIGEN_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "DRB345_SUBLOCUS_MISMATCH_AND_ANTIGEN_MISMATCH"
        )

    def test_LOCUS_MISMATCH_and_ALLELE_MISMATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.ALLELE_MISMATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(
            result, "DRB345_SUBLOCUS_MISMATCH_AND_ALLELE_MISMATCH"
        )

    def test_LOCUS_MISMATCH_and_ARD_MATCH(self):
        level1 = AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH
        level2 = AlleleMatchLevel.ARD_MATCH
        result = dummy_MatchResult._calculate_loci_match_high_resolution(
            level1, level2
        )
        self.assertEqual(result, "PARTIAL_DRB345_SUBLOCUS_MISMATCH")

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
                AlleleMatchLevel.ARD_MATCH,
                AlleleMatchLevel.ARD_MATCH
            )
        )
        self.assertEqual(
            result[1].allele_match_levels,
            (
                AlleleMatchLevel.ARD_MATCH,
                AlleleMatchLevel.ARD_MATCH
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
                "matching will be reported as NOT_ASSESSABLE",
                cm.output[0],
            )

        self.assertEqual(len(result), 2)
        self.assertEqual(
            result[0].allele_match_levels,
            (
                AlleleMatchLevel.ARD_MATCH,
                AlleleMatchLevel.ARD_MATCH
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
                AlleleMatchLevel.NOT_ASSESSABLE,
                AlleleMatchLevel.NOT_ASSESSABLE
            ),
        )


if __name__ == "__main__":
    unittest.main()
