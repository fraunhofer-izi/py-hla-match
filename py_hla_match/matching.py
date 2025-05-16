from py_hla_match.external import is_permissive_dpb1_match
import logging
from py_hla_match.models import HLAPair, Individual

from py_hla_match.hla import HLA
from enum import IntEnum
from typing import List, Tuple

from py_hla_match.exceptions import InvalidLocusComparisonError

logger = logging.getLogger(__name__)


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
            allele_match_levels: Tuple[AlleleMatchLevel, AlleleMatchLevel],
            compute_dpb1_permissive: bool = False
    ) -> None:

        self.patient = patient
        self.donor = donor
        self.allele_score = score
        self.allele_match_levels = allele_match_levels
        self.is_homozygous_patient = (
                self.patient.hla1.ard_redux_allele_string
                == self.patient.hla2.ard_redux_allele_string
        )
        self.dpb1_permissive = None
        if compute_dpb1_permissive:
            self._compute_dpb1_permissive()


    def get_match_level_for_resolution(self, resolution: str) -> str:
        """
        Get match level for a given resolution

        Args:
            resolution (str): Resolution level (basic, high, full)

        Returns:
            str: Match level for the given resolution
        """
        if resolution == "basic":
            return self.loci_match_basic_resolution
        elif resolution == "high":
            return self.loci_match_high_resolution
        elif resolution == "full":
            return self.loci_match_full_resolution
        else:
            raise ValueError(
                f"Unknown resolution level: {resolution}\n"
                f"Expected 'basic', 'high', or 'full'."
            )


    @property
    def loci_match_basic_resolution(self):
        if not hasattr(self, '_locus_match_basic_resolution'):
            self._locus_match_basic_resolution = \
                self._loci_level_match('basic_resolution')
        return self._locus_match_basic_resolution

    @property
    def loci_match_high_resolution(self):
        if not hasattr(self, '_locus_match_high_resolution'):
            self._locus_match_high_resolution = \
                self._loci_level_match('high_resolution')
        return self._locus_match_high_resolution

    @property
    def loci_match_full_resolution(self):
        if not hasattr(self, '_locus_match_full_resolution'):
            self._locus_match_full_resolution = \
                self._loci_level_match('full_resolution')
        return self._locus_match_full_resolution

    def _get_details(self) -> str:
        """
        TODO: not implemented yet
        """
        raise NotImplementedError("Not implemented yet.")

    def _loci_level_match(self, resolution):
        """
        TODO: base on clinician feedback
        """
        match_level_1, match_level_2 = self.allele_match_levels
        if resolution == 'basic_resolution':
            return self._calculate_loci_match_basic_resolution(
                match_level_1, match_level_2
            )

        elif resolution == 'high_resolution':
            return self._calculate_loci_match_high_resolution(
                match_level_1, match_level_2
            )

        elif resolution == 'full_resolution':
            raise NotImplementedError(
                "\'full_resolution\' not implemented yet.\n"
                "Please use either \'basic_resolution\' or "
                "\'high_resolution\'."
            )

        else:
            raise ValueError(
                f"Unknown resolution level: {resolution}\n"
                f"Expected 'basic_resolution', 'high_resolution', or "
                f"'full_resolution'."
            )

    def _calculate_loci_match_basic_resolution(
            self, match_level_1, match_level_2
    ):
        """
        TODO: base on clinician feedback
        Determines the basic resolution match status based on the allele match
        levels.

        Returns:
            str: "ARD_MATCH", "PARTIAL_ARD_MISMATCH", or "ARD_MISMATCH"
        """
        # type check
        if not all(
                isinstance(level, AlleleMatchLevel) for
                level in [match_level_1, match_level_2]
        ):
            raise TypeError(
                f"match_level_1 and match_level_2 must be instances of "
                f"{AlleleMatchLevel}, not {type(match_level_1)} and "
                f"{type(match_level_2)}."
            )
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

    def _calculate_loci_match_high_resolution(
            self, match_level_1, match_level_2
    ):
        """
        TODO: base on clinician feedback
        Determines the high resolution match status with detailed mismatch
        types.

        Returns:
            str: A string indicating the match status with high resolution
            mismatches.
        """
        # type check
        if not all(
                isinstance(level, AlleleMatchLevel) for
                level in [match_level_1, match_level_2]
        ):
            raise TypeError(
                f"match_level_1 and match_level_2 must be instances of "
                f"{AlleleMatchLevel}, not {type(match_level_1)} and "
                f"{type(match_level_2)}."
            )

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
            # resolve high resolution mismatch level
            return f"PARTIAL_{match_level_2.name}"
        elif (
                match_level_1 in mismatch_levels and
                match_level_2 in match_levels
        ):
            # resolve high resolution mismatch level
            return f"PARTIAL_{match_level_1.name}"
        # Both alleles are high resolution mismatch level
        elif (
                match_level_1 in mismatch_levels and
                match_level_2 in mismatch_levels and
                match_level_1 < match_level_2  # Order of mismatch "severity"
        ):
            return f"{match_level_1.name}_AND_{match_level_2.name}"
        elif (
                match_level_1 in mismatch_levels and
                match_level_2 in mismatch_levels and
                match_level_1 > match_level_2  # Order of mismatch "severity"
        ):
            return f"{match_level_2.name}_AND_{match_level_1.name}"
        # Additional sanity check
        elif (
                match_level_1 in mismatch_levels and
                match_level_2 in mismatch_levels and
                match_level_1 is match_level_2
        ):
            return f"DOUBLE_{match_level_1.name}"  # TODO: discuss terminology
        else:
            raise ValueError(
                f"Unexpected match levels {match_level_1.name}"
                f"and {match_level_2.name}"
            )

    def _compute_dpb1_permissive(self) -> None:
        # only compute for DPB1 locus
        if self.patient.locus != "DPB1":
            pass
        else:
            patient_dpb1 = self.patient.hla1.ard_redux_allele_string
            patient_dpb2 = self.patient.hla2.ard_redux_allele_string
            donor_dpb1 = self.donor.hla1.ard_redux_allele_string
            donor_dpb2 = self.donor.hla2.ard_redux_allele_string
            dpb1_permissive = is_permissive_dpb1_match(patient_dpb1, patient_dpb2, donor_dpb1, donor_dpb2)
            self.dpb1_permissive = dpb1_permissive


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


def allele_pair_match(patient: HLAPair, donor: HLAPair, calculate_dpb1_permissive: bool = False) -> MatchResult:
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
        calculate_dpb1_permissive=calculate_dpb1_permissive
    )



def multi_locus_match(
        patient: Individual,
        donor: Individual
) -> List[MatchResult]:
    """
    Calculate compatibility of patient and donor for all recorded patient loci

    Args:
        patient (Individual): Patient object
        donor (Individual): Donor object

    Returns:
        List[MatchResult]: List of MatchResult objects for each locus
    """
    results = []

    # create dict
    donor_hla_dict = {hla_pair.locus: hla_pair for hla_pair in donor.hla_data}

    for patient_hla_pair in patient.hla_data:
        # check if there is a matching locus in the donor data
        locus = patient_hla_pair.locus

        if locus in donor_hla_dict:
            results.append(
                allele_pair_match(patient_hla_pair, donor_hla_dict[locus])
            )
        else:
            logger.warning(
                f"Locus {locus} not found in donor data and will be excluded"
                f" from the results."
            )

    return results
