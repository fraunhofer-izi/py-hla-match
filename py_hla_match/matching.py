from py_hla_match.models import HLAPair
from py_hla_match.hla import HLA
from enum import IntEnum
from typing import Tuple

from py_hla_match.exceptions import InvalidLocusComparisonError


class AlleleMatchLevel(IntEnum):
    """
    Following hla nomenclature:
    LOCUS_MISMATCH: Mismatch of a particular HLA locus i.e. DRB1
    ALLELE_GROUP_MISMATCH:
    ALLELE_MISMATCH = 2
    ARD_MATCH = 3
    SYNONYMOUS_VARIANT_MATCH = 4
    NON_CODING_VARIANT_MATCH = 5    
    cf.https://hla.alleles.org/nomenclature/naming.html
    """
    LOCUS_MISMATCH = 0
    ALLELE_GROUP_MISMATCH = 1
    ALLELE_MISMATCH = 2
    ARD_MATCH = 3
    SYNONYMOUS_VARIANT_MATCH = 4
    NON_CODING_VARIANT_MATCH = 5


class MatchResult:
    """
    Report class on patient/donor compatibility for a single HLA locus

    Attributes:
        patient (HLAPair): Patient object (HLA alleles)
        donor (HLAPair): Donor object (HLA alleles)
        score (int): (internal) score of correct allele pairing (patient/donor)
        allele_matches (Tuple[AlleleMatchLevel, AlleleMatchLevel]): Match
        results
    """

    def __init__(
            self,
            patient: HLAPair,
            donor: HLAPair,
            score: int,
            allele_match_levels: Tuple[AlleleMatchLevel, AlleleMatchLevel]
    ) -> None:

        self.patient = patient
        self.donor = donor
        self.allele_score = score
        self.allele_match_levels = allele_match_levels
        self.is_homozygous_patient = (
                   self.patient.hla1.ard_redux_allele_string
                   == self.patient.hla2.ard_redux_allele_string
               )

    @staticmethod
    def _get_details() -> str:
        """
        TODO: not implemented yet
        """
        return None

    @staticmethod
    def _loci_level_match():
        """
        TODO: base on clinician feedback
        """
        return None


def allele_match(hla1: HLA, hla2: HLA) -> AlleleMatchLevel:
    """
    Compares two HLA alleles and returns a MatchLevel

    Args:
        hla1: First HLA allele object
        hla2: Second HLA allele object

    Returns:
        MatchLevel enum value indicating position of matches and mismatch
        (cf. HLA nomenclature)
    Raises:
        TypeError: If hla1 or hla2 is not an instance of HLA
        InvalidLocusComparisonError: If hla1 and hla2 have incompatible loci
    """

    if not isinstance(hla1, HLA):
        raise TypeError(
            f"hla1 must be an instance of HLA, not {type(hla1).__name__}."
        )
    if not isinstance(hla2, HLA):
        raise TypeError(
            f"hla2 must be an instance of HLA, not {type(hla2).__name__}."
        )

    if hla1.locus != hla2.locus:
        if (
            'DRB' in hla1.locus and 'DRB' in hla2.locus and
            'DRB1' not in hla1.locus and 'DRB1' not in hla2.locus
        ):
            return AlleleMatchLevel.LOCUS_MISMATCH
        else:
            raise InvalidLocusComparisonError(hla1.locus, hla2.locus)

    if hla1.ard_redux_allele_group != hla2.ard_redux_allele_group:
        return AlleleMatchLevel.ALLELE_GROUP_MISMATCH

    if hla1.ard_redux_allele != hla2.ard_redux_allele:
        return AlleleMatchLevel.ALLELE_MISMATCH

    # from here on we have at least an ARD level match

    # Check for group code
    if (
            hla1.group_code is not None
            or hla2.group_code is not None
    ):
        # If group_code is provided, we will not exceed ARD level match
        return AlleleMatchLevel.ARD_MATCH

    # Check for suffix
    if (
            hla1.suffix is not None
            or
            hla2.suffix is not None
    ):
        # TODO: implement logic for suffixes
        # for now, skip these alleles
        return AlleleMatchLevel.ARD_MATCH

    # Compare specific allele
    if hla1.allele != hla2.allele:
        return AlleleMatchLevel.ARD_MATCH  # ARD level match remains

    # To save some time continue with synonymous variant+ only if available
    if (
        hla1.synonymous_variant is not None
        and
        hla2.synonymous_variant is not None
    ):
        if hla1.synonymous_variant == hla2.synonymous_variant:
            # continue with non coding variant
            if (
                    hla1.non_coding_variant is not None
                    and
                    hla2.non_coding_variant is not None
            ):
                if hla1.non_coding_variant == hla2.non_coding_variant:
                    # highest resulution match
                    return AlleleMatchLevel.NON_CODING_VARIANT_MATCH
                else:
                    # synonymous variant match remains
                    return AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
            else:
                # synonymous variant match remains
                return AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        else:
            # ARD level match remains, if synonymous variant(s) are not equal
            return AlleleMatchLevel.ARD_MATCH
    else:
        # If no information on synonymous variant, ARD level match remains
        return AlleleMatchLevel.ARD_MATCH


def _get_correct_allele_pairing(
    patient_alleles: HLAPair, donor_alleles: HLAPair
) -> Tuple[int, Tuple[AlleleMatchLevel, AlleleMatchLevel]]:
    """
    Determines the correct pairing of patient and donor HLA allele by
    evaluating all possible combinations

    Args:
        patient_alleles: Tuple[HLA, HLA]: Tuple of two patient HLA alleles
        donor_alleles: Tuple[HLA, HLA]: Tuple of two donor HLA alleles

    Returns:
        Tuple[int, Tuple[AlleleMatchLevel, AlleleMatchLevel]]:
            - best_score (int): Best score among all possible allele pairings
            - correct_pairing (Tuple[AlleleMatchLevel, AlleleMatchLevel]):
            Tuple containing match levels of correct allele pairing

    Notes:
        - The function assumes that both `patient_alleles` and `donor_alleles`
        contain exactly two alleles each
        - Considers two possible pairings:
            1. (patient_hla1, donor_hla1) and (patient_hla2, donor_hla2)
            2. (patient_hla1, donor_hla2) and (patient_hla2, donor_hla1)
        - Best score-pairing is returned
        - If pairings return equal score, the first pairing is returned
    """
    pairings = [
        (
            patient_alleles.hla1, donor_alleles.hla1,
            patient_alleles.hla2, donor_alleles.hla2
        ),
        (
            patient_alleles.hla1, donor_alleles.hla2,
            patient_alleles.hla2, donor_alleles.hla1
        ),
    ]

    best_score = -1
    correct_pairing = None

    for pairing in pairings:
        patient_hla1, donor_hla1, patient_hla2, donor_hla2 = pairing

        # Compute match levels
        allele_match1 = allele_match(patient_hla1, donor_hla1)
        allele_match2 = allele_match(patient_hla2, donor_hla2)

        # Calculate total score
        score = allele_match1 + allele_match2

        # Update best score and (if better) pairing
        if score > best_score:
            best_score = score
            correct_pairing = [allele_match1, allele_match2]

    return best_score, correct_pairing


def allele_pair_match(patient: HLAPair, donor: HLAPair) -> MatchResult:
    """
    Matching of two patient and donor HLA alleles encoding the HLA gene

    Args:
        patient (Patient): Patient object containing two HLA alleles
        donor (Donor): Donor object containing two HLA alleles

    Returns:
        MatchResult: Class storing matching results

    Notes:
        - The function assumes that both, patient and donor, have exactly two
        HLA alleles
        - Uses `get_correct_allele_pairing` function to evaluate all possible
        allele pairings and selects the one with the highest score
    """
    # Get correct allele pairing and its score
    score, correct_pairing = _get_correct_allele_pairing(patient, donor)

    # Return match result
    return MatchResult(
        patient=patient,
        donor=donor,
        score=score,
        allele_match_levels=correct_pairing,
    )
