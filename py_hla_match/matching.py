from py_hla_match.models import Patient, Donor
from py_hla_match.hla import HLA
from enum import IntEnum
from typing import Tuple

from py_hla_match.exceptions import InvalidLocusComparisonError


class AlleleMatchLevel(IntEnum):
    """
    Following hla nomenclature:
    LOCUS_MISMATCH: Mismatch at a particular HLA locus
    ALLELE_GROUP_MISMATCH: Mismatch at the goup code
    ALLELE_MISMATCH: Mismatch at the allele level
    ARD_MATCH: ARD level match
    SYNONYMOUS_VARIANT_MATCH: Synonymous variant match
    NON_CODING_VARIANT_MATCH: Non-coding variant match
    cf.https://hla.alleles.org/nomenclature/naming.html
    """
    LOCUS_MISMATCH = 0
    ALLELE_GROUP_MISMATCH = 1
    ALLELE_MISMATCH = 2
    ARD_MATCH = 3
    SYNONYMOUS_VARIANT_MATCH = 4
    NON_CODING_VARIANT_MATCH = 5


class MatchResult():
    """
    Report class on patient/donor compatibility for an HLA locus

    Attributes:
        patient (Patient): Patient object (HLA alleles)
        donor (Donor): Donor object (HLA alleles)
        score (int): (internal) score of correct allele pairing (patient/donor)
        allele_matches (Tuple[AlleleMatchLevel, AlleleMatchLevel]): Match
        results
    """

    def __init__(
            self,
            patient: Patient,
            donor: Donor,
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

    @property
    def loci_match_basic_resolution(self):
        if not hasattr(self, '_locus_match_basic_resolution'):
            self._locus_match_basic_resolution = \
                self._calculate_locus_match('basic_resolution')
        return self._locus_match_basic_resolution

    @property
    def loci_match_high_resolution(self):
        if not hasattr(self, '_locus_match_high_resolution'):
            self._locus_match_high_resolution = \
                self._calculate_locus_match('high_resolution')
        return self._locus_match_high_resolution

    @property
    def loci_match_full_resolution(self):
        if not hasattr(self, '_locus_match_high_resolution'):
            self._locus_match_high_resolution = \
                self._calculate_locus_match('high_resolution')
        return self._locus_match_high_resolution

    def _get_details(self) -> str:
        """
        TODO: not implemented yet
        """
        return None

    def _loci_level_match(self, resolution):
        """
        TODO: base on clinician feedback
        """
        match_level_1, match_level_2 = self.allele_match_levels
        if resolution == 'basic_resolution':
            return self._calculate_basic_resolution(match_level_1, match_level_2)

        elif resolution == 'high_resolution':
            return self._calculate_high_resolution(match_level_1, match_level_2)

        else:
            raise ValueError(f"Unknown resolution level: {resolution}")

    def _calculate_loci_match_basic_resolution(self, match_level_1, match_level_2):
        """
        Determines the basic resolution match status based on the allele match
        levels.

        Returns:
            str: "ARD_MATCH", "PARTIAL_ARD_MISMATCH", or "ARD_MISMATCH"
        """
        # Group AlleleMatchLevels into basic resolution match and mismatch
        # levels
        match_levels = {
            AlleleMatchLevel.ARD_MATCH,
            AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH,
            AlleleMatchLevel.NON_CODING_VARIANT_MATCH,
        }
        mismatch_levels = {
            AlleleMatchLevel.LOCUS_MISMATCH,
            AlleleMatchLevel.ALLELE_GROUP_MISMATCH,
            AlleleMatchLevel.ALLELE_MISMATCH,
        }

        # Combine grouped basic resolution match level to "MATCH"
        if match_level_1 in match_levels and match_level_2 in match_levels:
            return "ARD_MATCH"
        # Partial mismatch if one allele is basic resolution match level and
        # the other is basic resolution mismatch level
        elif (
            match_level_1 in match_levels and
            match_level_2 in mismatch_levels
        ) or (
            match_level_1 in mismatch_levels and
            match_level_2 in match_levels
        ):
            return "PARTIAL_ARD_MISMATCH"
        # Both alleles are basic mismatch level
        else:
            return "ARD_MISMATCH"

    def _calculate_loci_match_high_resolution(self, match_level_1, match_level_2):
        """
        Determines the high resolution match status with detailed mismatch
        types.

        Returns:
            str: A string indicating the match status with high resolution
            mismatches.
        """

        # Group AlleleMatchLevels into high resolution match and mismatch
        # levels
        match_levels = {
            AlleleMatchLevel.ARD_MATCH,
            AlleleMatchLevel.SYNONYMOUS_VARIANT_MATCH,
            AlleleMatchLevel.NON_CODING_VARIANT_MATCH,
        }
        mismatch_levels = {
            AlleleMatchLevel.LOCUS_MISMATCH,
            AlleleMatchLevel.ALLELE_GROUP_MISMATCH,
            AlleleMatchLevel.ALLELE_MISMATCH,
        }

        # Combine grouped high resolution match level to "MATCH"
        if match_level_1 in match_levels and match_level_2 in match_levels:
            return "ARD_MATCH"
        # Partial mismatch if one allele is high resolution match level and
        # the other is high resolution mismatch level
        elif (
            match_level_1 in match_levels and
            match_level_2 in mismatch_levels
        ):
            return f"PARTIAL_{match_level_2}"
        elif (
            match_level_1 in mismatch_levels and
            match_level_2 in match_levels
        ):
            return f"PARTIAL_{match_level_1}"
        # Both alleles are high resolution mismatch level
        elif (
            match_level_1 in mismatch_levels and
            match_level_2 in mismatch_levels and
            match_level_1 < match_level_2  # Order of mismatch "severity"
        ):
            return f"{match_level_1}_AND_{match_level_2}"
        elif (
            match_level_1 in mismatch_levels and
            match_level_2 in mismatch_levels and
            match_level_1 > match_level_2  # Order of mismatch "severity"
        ):
            return f"{match_level_2}_AND_{match_level_1}"
        # Additional sanity check
        elif (
            match_level_1 in mismatch_levels and
            match_level_2 in mismatch_levels and
            match_level_1 is match_level_2
        ):
            return f"DOUBLE_{match_level_1}"  # TODO: discuss terminology
        else:
            raise ValueError(
                f"Unexpected match levels {match_level_1} and {match_level_2}"\
                f"in {__name__}"
            )


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
    patient_alleles: Tuple[HLA, HLA], donor_alleles: Tuple[HLA, HLA]
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
    patient_alleles = [patient.hla1, patient.hla2]
    donor_alleles = [donor.hla1, donor.hla2]

    # Get correct allele pairing and its score
    score, correct_pairing = _get_correct_allele_pairing(
        patient_alleles, donor_alleles
    )

    # Return match result
    return MatchResult(
        patient=patient,
        donor=donor,
        score=score,
        allele_match_levels=correct_pairing,
    )
