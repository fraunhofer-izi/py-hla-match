import logging
from py_hla_match.models import HLAPair, Individual
from py_hla_match.hla import HLA
from typing import List, Tuple, Optional

from py_hla_match.policy import AlleleMatchLevel, ExpressionSuffixMatchLevel
from py_hla_match.config import (
    get_config,
    HLAMatchConfig,
)
from py_hla_match.exceptions import InvalidLocusComparisonError
from py_hla_match.external import DPB1TCEStatus, query_dpb1_tce

logger = logging.getLogger(__name__)


class MatchResult:
    """
    Result object for comparing two HLA genotype pairs at a single locus.

    This class is designed for research use to describe HLA match or
    mismatch categories between two individuals.

    :ivar patient: HLA allele pair in the 'patient' role
    :ivar donor: HLA allele pair in the 'donor' role
    :ivar pairing_score: internal, ordinal score summarising the two
        AlleleMatchLevel values.
    """

    def __init__(
            self,
            patient: HLAPair,
            donor: HLAPair,
            pairing_score: int,
            allele_match_levels: Tuple[AlleleMatchLevel, AlleleMatchLevel],
    ) -> None:

        self.patient = patient
        self.donor = donor
        self.pairing_score = pairing_score
        self.allele_match_levels = allele_match_levels

        # optional external matching information
        self.dpb1_tce_status: Optional[DPB1TCEStatus] = None

        # check homozygous patient
        self.is_homozygous_patient = (
            # get boolean if patient alleles are equal
            (
                self.patient.hla1.ard_redux_allele_string
                == self.patient.hla2.ard_redux_allele_string
            )
            # if ard redux is available
            if (self.patient.hla1.ard_redux_allele_string
                and self.patient.hla2.ard_redux_allele_string)
            # else None
            else None
        )

    def get_match_level_for_resolution(self, resolution: str) -> str:
        """
        Get match level for a given resolution

        :param resolution: Resolution level (basic, high, full)
        :return: Match level for the given resolution
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
        TODO: refine resolution categories based on domain-expert input
        """
        match_level_1, match_level_2 = self.allele_match_levels

        if (
            match_level_1 is AlleleMatchLevel.NOT_APPLICABLE
            or match_level_2 is AlleleMatchLevel.NOT_APPLICABLE
        ):
            return AlleleMatchLevel.NOT_APPLICABLE.name

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
        TODO: base on domain expert
        Determines the basic resolution match status based on the allele match
        levels.

        :return: "ARD_MATCH", "PARTIAL_ARD_MISMATCH", or "ARD_MISMATCH"
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
        TODO: base on domain expert
        Determines the high resolution match status with detailed mismatch
        types.

        :return: A string indicating the match status with high resolution
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

    @staticmethod
    def _api_allele(hla: HLA) -> Optional[str]:
        """
        Get highest resolution string to query EBI TCE API.
        """
        if hla.ard_redux_allele_string:
            return hla.ard_redux_allele_string
        if hla.allele_group:
            return f"{hla.locus}*{hla.allele_group}"
        return None

    def get_dpb1_tce_status(
            self,
            api_version: str = "3.0",
            timeout: int = 10
    ) -> Optional[DPB1TCEStatus]:
        """
        Calculate DPB1 permissive/non-permissive classification via EBI API.

        Intended for research workflows.

        WARNING: may slow things down significantly.

        Sets self.dpb1_permissive to one of:
        - DPB1TCEStatus
        """
        if self.patient.locus != "DPB1":
            logger.debug(
                f"Not applicable to {self.patient.locus}."
            )
            return None

        patient_dpb1 = self._api_allele(self.patient.hla1)
        patient_dpb2 = self._api_allele(self.patient.hla2)
        donor_dpb1 = self._api_allele(self.donor.hla1)
        donor_dpb2 = self._api_allele(self.donor.hla2)

        if not all(
            [patient_dpb1, patient_dpb2, donor_dpb1, donor_dpb2]
        ):
            logger.warning(
                f"One or more required alleles are missing for DPB1 to call "
                f"EBI API, got P1:'{patient_dpb1}', P2:'{patient_dpb2}', "
                f"D1:'{donor_dpb1}', D2:'{donor_dpb2}'. dpb1_tce_status "
                f"remains unchanged ({self.dpb1_tce_status})"
            )
            return None

        # if we are here, the query should be valid
        dpb1_tce_status = query_dpb1_tce(
            patient_dpb1=patient_dpb1,
            patient_dpb2=patient_dpb2,
            donor_dpb1=donor_dpb1,
            donor_dpb2=donor_dpb2,
            version=api_version,
            timeout=timeout
        )
        # update match result and return DPB1TCEStatus
        self.dpb1_tce_status = dpb1_tce_status
        return self.dpb1_tce_status


def _map_expression_decision(
    decision: ExpressionSuffixMatchLevel,
) -> Optional[AlleleMatchLevel]:
    if decision is ExpressionSuffixMatchLevel.IGNORE:
        return None
    mapping = {
        ExpressionSuffixMatchLevel.NOT_APPLICABLE:
            AlleleMatchLevel.NOT_APPLICABLE,
        ExpressionSuffixMatchLevel.ALLELE_MISMATCH:
            AlleleMatchLevel.ALLELE_MISMATCH,
        ExpressionSuffixMatchLevel.ALLELE_GROUP_MISMATCH:
            AlleleMatchLevel.ALLELE_GROUP_MISMATCH,
        ExpressionSuffixMatchLevel.ARD_MATCH:
            AlleleMatchLevel.ARD_MATCH,
    }
    return mapping[decision]


def _apply_expression_suffix_policy(
    hla1: HLA, hla2: HLA, cfg: HLAMatchConfig
) -> Optional[AlleleMatchLevel]:
    """
    Apply configurable expression-suffix policy once ARD is equivalent.
    """
    suffix1, suffix2 = hla1.suffix, hla2.suffix
    if suffix1 is None and suffix2 is None:
        return None
    rules = cfg.expression_suffix_policy
    # Any 'Q' present (defaults to NOT_APPLICABLE)
    if (
            (suffix1 in rules.ambiguous_suffixes) or
            (suffix2 in rules.ambiguous_suffixes)
    ):
        return _map_expression_decision(rules.q_present)
    # Any risk suffixes present
    risk = rules.risk_suffixes
    risk1 = suffix1 in risk if suffix1 is not None else False
    risk2 = suffix2 in risk if suffix2 is not None else False
    if risk1 and risk2:
        if suffix1 == suffix2:
            return _map_expression_decision(rules.equal_risk)
        return _map_expression_decision(rules.risk_vs_different_risk)
    if (risk1 and suffix2 is None) or (risk2 and suffix1 is None):
        return _map_expression_decision(rules.risk_vs_none)
    return None


def allele_match(hla1: HLA, hla2: HLA) -> AlleleMatchLevel:
    """
    Compares two HLA alleles and returns a MatchLevel

    :param hla1: First HLA allele object
    :param hla2: Second HLA allele object

    :return: MatchLevel enum value indicating position of matches and mismatch
        (cf. HLA nomenclature)
    :raises TypeError: If hla1 or hla2 is not an instance of HLA
    :raises InvalidLocusComparisonError: If hla1 and hla2 have incompatible loci
    """

    if not isinstance(hla1, HLA):
        raise TypeError(
            f"hla1 must be an instance of HLA, not {type(hla1).__name__}."
        )
    if not isinstance(hla2, HLA):
        raise TypeError(
            f"hla2 must be an instance of HLA, not {type(hla2).__name__}."
        )

    # (1) ARD COMPARISON

    # first check if loci match (NOTE: DRB3/4/5 hard coded to locus DRB345)
    if hla1.locus != hla2.locus:
        raise InvalidLocusComparisonError(hla1.locus, hla2.locus)

    # for locus DRB345, we still stored the original DRB3/4/5 sub-locus
    if hla1.locus == 'DRB345' and hla1.drb_sub_locus != hla2.drb_sub_locus:
        return AlleleMatchLevel.LOCUS_MISMATCH

    if min(
        hla1.has_resolution_level(), hla2.has_resolution_level()
    ) < 1:
        # no allele fields
        return AlleleMatchLevel.NOT_APPLICABLE

    if min(
        hla1.has_resolution_level(), hla2.has_resolution_level()
    ) < 2:
        # check if allele groups differ
        if hla1.allele_group != hla2.allele_group:
            return AlleleMatchLevel.ALLELE_GROUP_MISMATCH
        # else we cannot determine a match level (missing data)
        else:
            return AlleleMatchLevel.NOT_APPLICABLE

    # --- from here on we have at least two-field resolution ---
    if (
        hla1.ard_redux_allele_string is None
        or hla2.ard_redux_allele_string is None
    ):
        # NOTE: this should never happen (!)
        raise RuntimeError(
            f"HLA parsing failed for '{hla1.allele_string}' or "
            f"'{hla2.allele_string}'. Please report this issue."
        )

    if hla1.ard_redux_allele_group != hla2.ard_redux_allele_group:
        # NOTE: this should never happen (!) since we check allele_group above
        return AlleleMatchLevel.ALLELE_GROUP_MISMATCH

    if hla1.ard_redux_allele != hla2.ard_redux_allele:
        return AlleleMatchLevel.ALLELE_MISMATCH

    # (2) EXPRESSION COMPARISON (suffixes)

    # Check for suffix
    if (
            hla1.suffix is not None
            or hla2.suffix is not None
    ):
        suffix_level = _apply_expression_suffix_policy(
            hla1, hla2, get_config()
        )
        if suffix_level is not None:
            return suffix_level

    # from here on we have at least an ARD level match that is NOT effected by
    # expression differences (suffixes)

    # (3) FULL ALLELE COMPARISON

    # Check for group code
    if (
            hla1.group_code is not None
            or hla2.group_code is not None
    ):
        # If P group_code is provided, we will not exceed ARD level match
        return AlleleMatchLevel.ARD_MATCH
        # If G group_code is provided, we could exceed ARD level match
        # TODO: continue

    # Compare specific allele (we check this again to continue with
    # synonymous variant+ only if non-redux allele strings are equal)
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

    :param patient_alleles: Tuple[HLA, HLA]: Tuple of two patient HLA alleles
    :param donor_alleles: Tuple[HLA, HLA]: Tuple of two donor HLA alleles

    :return: Tuple[int, Tuple[AlleleMatchLevel, AlleleMatchLevel]]:
        - best_score (int): Best pairing score of correct allele pairing
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

    best_score = float('-inf')  # lowest possible score, cf. AlleleMatchLevel
    correct_pairing = None

    for pairing in pairings:
        patient_hla1, donor_hla1, patient_hla2, donor_hla2 = pairing

        # Compute match levels
        allele_match1 = allele_match(patient_hla1, donor_hla1)
        allele_match2 = allele_match(patient_hla2, donor_hla2)

        # Calculate pairing score
        pairing_score = allele_match1 + allele_match2

        # Update best score and (if better) pairing
        if pairing_score > best_score:
            best_score = pairing_score
            correct_pairing = (allele_match1, allele_match2)

    return best_score, correct_pairing


def allele_pair_match(patient: HLAPair, donor: HLAPair) -> MatchResult:
    """
    Compute research match/mismatch levels for two HLA allele pairs, one in
    the 'patient' role and one in the 'donor' role.

    Intended for research workflows.

    :param patient: Patient object containing two HLA alleles
    :param donor: Donor object containing two HLA alleles

    :return: Class storing matching results

    Notes:
        
                - The function assumes that both, patient and donor, have exactly two
                    HLA alleles
                - Uses `_get_correct_allele_pairing` function to evaluate all possible
                    allele pairings and selects the one with the highest score (i.e.
                    correct pairing)
    """
    # Get correct allele pairing and its score
    pairing_score, correct_pairing = _get_correct_allele_pairing(
        patient, donor
    )

    # Create MatchResult object
    match_result = MatchResult(
        patient=patient,
        donor=donor,
        pairing_score=pairing_score,
        allele_match_levels=correct_pairing
    )

    # Return match result
    return match_result


def multi_locus_match(
        patient: Individual,
        donor: Individual
) -> List[MatchResult]:
    """
    Compute HLA match/mismatch categories between two Individuals for all loci
    that are typed in the first Individual.

    Intended for research workflows.

    :param patient: Patient object
    :param donor: Donor object

    :return: List of MatchResult objects for each locus
    """
    results: List[MatchResult] = []

    # quick lookup for donor pairs by locus
    donor_dict = {pair.locus: pair for pair in donor.hla_data}

    for patient_pair in patient.hla_data:
        locus = patient_pair.locus

        if locus in donor_dict:
            donor_pair = donor_dict[locus]
        else:
            logger.warning(
                f"Locus {locus} not found in donor data – "
                "matching will be reported as NOT_APPLICABLE."
            )
            if locus == "DRB345":
                donor_pair = HLAPair(HLA("DRBX*NA"), HLA("DRBX*NA"))
            else:
                donor_pair = HLAPair(HLA(f"{locus}*NA"), HLA(f"{locus}*NA"))

        # compute match (missing donor pair will propagate NOT_APPLICABLE)
        match_result = allele_pair_match(patient_pair, donor_pair)

        # additional diagnostics
        if all(level == AlleleMatchLevel.NOT_APPLICABLE
               for level in match_result.allele_match_levels):
            logger.warning(
                f"Typing resolution insufficient for locus {locus} "
                f"(patient {patient_pair} / donor {donor_pair})."
            )
        results.append(match_result)

    return results
