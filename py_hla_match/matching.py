from py_hla_match.models import Patient, Donor
from py_hla_match.hla import HLA
from enum import IntEnum
from typing import List, Tuple

from py_hla_match.exceptions import InvalidLocusComparisonError


class AlleleMatchLevel(IntEnum):
    LOCUS_MISMATCH = 0
    ALLELE_GROUP_MISMATCH = 1
    ALLELE_MISMATCH = 2
    ARD_MATCH = 3
    SYNONYMOUS_VARIANT_MATCH = 4
    NON_CODING_VARIANT_MATCH = 5


class MatchResult():
    """
    Report class on patient/donor compatibility.

    Attributes:
        patient (Patient): Patient object (HLA alleles)
        donor (Donor): Donor object (HLA alleles)
        score (int): Score of correct allele pairing (patient/donor).
        allele_matches (List[AlleleMatchLevel]): List of match levels for the
        correct allele pairing.
    """

    def __init__(
            self,
            patient: Patient,
            donor: Donor,
            score: int,
            allele_match_levels: List[AlleleMatchLevel]
    ) -> None:

        self.patient = patient
        self.donor = donor
        self.allele_score = score
        self.allele_match_levels = allele_match_levels
        self.is_homozygous_patient = (
                   self.patient.hla1.ard_redux_allele_string ==
                   self.patient.hla2.ard_redux_allele_string
               )

    def _get_details(self) -> str:
        """
        TODO: not implemented yet
        """
        return None

    def _gene_level_match(self):
        """
        TODO: base on clinician feedback
        """
        return None


def allele_match(hla1, hla2) -> AlleleMatchLevel:
    """
    Compares two HLA alleles and returns a MatchLevel.

    Args:
        hla1: First HLA allele object.
        hla2: Second HLA allele object.

    Returns:
        MatchLevel enum value indicating position of matches and mismatch
        (cf. HLA nomenclature).
    """

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
    if hla1.group_code != hla2.group_code:
        return AlleleMatchLevel.ARD_MATCH  # TODO: handle 'G' groups

    # Check for suffix
    if hla1.suffix != hla2.suffix:
        return AlleleMatchLevel.ARD_MATCH  # TODO: handle suffixes

    # Compare specific allele
    if hla1.allele != hla2.allele:
        return AlleleMatchLevel.ARD_MATCH  # ARD level match remains

    # Continue with synonymous variant only if available
    if (
        hla1.synonymous_variant is not None
        and
        hla2.synonymous_variant is not None
    ):
        if hla1.synonymous_variant != hla2.synonymous_variant:
            return AlleleMatchLevel.ARD_MATCH
        else:
            return AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH

    # Continue with non coding variant only if available
    if (
            hla1.non_coding_variant is not None
            and
            hla2.non_coding_variant is not None
    ):
        if hla1.non_coding_variant != hla2.non_coding_variant:
            return AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH
        else:
            return AlleleMatchLevel.NON_CODING_VARIANT_MATCH

        # we arrived at highest resolution
        return AlleleMatchLevel.NON_CODING_VARIANT_MATCH
    else:
        return AlleleMatchLevel.ARD_MATCH  # ARD level match remains


def compare_allele_pairs(
        alleles: Tuple[HLA, HLA, HLA, HLA]
) -> Tuple[int, AlleleMatchLevel, AlleleMatchLevel]:
    """
    Compares pairs of HLA alleles between patient and donor and calculatess a
    score based on thier compatibility.

    Args:
        alleles (tuple): A tuple containing patient and donor HLA alleles in
        the following order:
            - alleles[0]: patient_hla1 (HLA)
            - alleles[1]: donor_hla1 (HLA)
            - alleles[2]: patient_hla2 (HLA)
            - alleles[3]: donor_hla2 (HLA)

    Returns:
        tuple:
            - total_score (int): The sum of the match levels for both allele
            pairs.
            - allele_match1 (AlleleMatchLevel): The match level for the first
            allele pair.
            - allele_match2 (AlleleMatchLevel): The match level for the second
            allele pair.
    """
    patient_hla1, donor_hla1, patient_hla2, donor_hla2 = alleles

    # matching
    allele_match1 = allele_match(patient_hla1, donor_hla1)
    allele_match2 = allele_match(patient_hla2, donor_hla2)

    # Calculate score
    total_score = allele_match1 + allele_match2

    return total_score, allele_match1, allele_match2


def get_correct_allele_pairing(
    patient_alleles: List[HLA], donor_alleles: List[HLA]
) -> Tuple[int, List[AlleleMatchLevel]]:
    """
    Determines the correct allele pairing between patient and a donor HLA
    alleles by evaluating all possible combinations and selecting the one with
    the highest compatibility score.

    Args:
        patient_alleles (List[HLA]): List of two patient HLA alleles
        donor_alleles (List[HLA]): List of two donor HLA alleles

    Returns:
        Tuple[int, List[AlleleMatchLevel]]:
            - best_score (int): The highest total score among all possible
            allele pairings
            - correct_pairing (List[AlleleMatchLevel]): List containing the
            match levels of the correct allele pairing.

    Notes:
        - The function assumes that both `patient_alleles` and `donor_alleles`
        contain exactly two alleles each.
        - Considers two possible pairings:
            1. (patient_hla1, donor_hla1) and (patient_hla2, donor_hla2)
            2. (patient_hla1, donor_hla2) and (patient_hla2, donor_hla1)
        - Pairing with the highest score is returned.
        - If both pairings return the same score, the first pairing is
        returned.
    """
    pairings = [
        (
            patient_alleles[0], donor_alleles[0],
            patient_alleles[1], donor_alleles[1]
        ),
        (
            patient_alleles[0], donor_alleles[1],
            patient_alleles[1], donor_alleles[0]
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


def allele_pair_match(patient: Patient, donor: Donor) -> MatchResult:
    """
    Determines the match result between patient and donor HLA gene.

    Args:
        patient (Patient): The patient object containing two HLA alleles.
        donor (Donor): The donor object containing two HLA alleles.

    Returns:
        MatchResult: Class storing the result of the gene matching function

    Notes:
        - The function assumes that both, patient and donor, have exactly two
        HLA alleles.
        - Uses `get_correct_allele_pairing` function to evaluate all possible
        allele pairings and selects the one with the highest score.
    """
    patient_alleles = [patient.hla1, patient.hla2]
    donor_alleles = [donor.hla1, donor.hla2]

    # Get correct allele pairing and its score
    score, correct_pairing = get_correct_allele_pairing(
        patient_alleles, donor_alleles
    )

    # Return match result
    return MatchResult(
        patient=patient,
        donor=donor,
        score=score,
        allele_match_levels=correct_pairing,
    )


if __name__ == "__main__":
    # Create HLA alleles for patient and donor
    patient_hla1 = HLA("B*07:02")
    patient_hla2 = HLA("B*44:02:01")
    patient = Patient(hla1=patient_hla1, hla2=patient_hla2)

    donor_hla1 = HLA("B*45:02:01G")
    donor_hla2 = HLA("B*18:01:01G")
    donor = Donor(hla1=donor_hla1, hla2=donor_hla2)

    # Perform gene matching
    match_result = allele_pair_match(patient, donor)

    # Access stored information
    print(f"Match Score: {match_result.allele_score}")
    print(f"Allele Match Levels: {match_result.allele_match_levels}")
    print(f"Is the patient homozygous? {match_result.is_homozygous_patient}")
