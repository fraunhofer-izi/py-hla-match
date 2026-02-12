import logging
from typing import List, Tuple, Optional
from dataclasses import dataclass

from py_hla_match.models import HLAPair, Individual
from py_hla_match.hla import HLA
from py_hla_match.policy import (
    AlleleMatchLevel,
    ExpressionSuffixMatchLevel,
    ARDMatchLevel,
    ARDMatchLevelCertainty,
    MolecularMatchLevel,
    MolecularMatchLevelCertainty
)
from py_hla_match.config import (
    get_config,
    HLAMatchConfig,
)
from py_hla_match.exceptions import (
    InvalidLocusComparisonError,
    ARDMatchRefinementError
)
from py_hla_match.external import DPB1TCEStatus, query_dpb1_tce
from py_hla_match.singleton import get_ard_instance


logger = logging.getLogger(__name__)


class MatchResult:
    """
    Result object for comparing two HLA genotype pairs at a single locus.

    This class is designed for research use to describe HLA match or
    mismatch categories between two individuals.

    :ivar patient: HLA allele pair in the 'patient' role.
    :ivar donor: HLA allele pair in the 'donor' role.
    :ivar pairing_score: Internal ordinal score summarising the two
        ``AlleleMatchLevel`` values.
    :ivar allele_match_levels: Tuple of ``AlleleMatchLevel`` values for the
        two allele-level comparisons (patient allele 1 vs donor allele X,
        patient allele 2 vs donor allele Y).
    :ivar ard_match_levels: Tuple of ``ARDMatchLevel`` values refining
        ARD-equivalent allele pairs (``NOT_APPLICABLE`` if not ARD-matched).
    :ivar ard_match_certainties: Tuple of ``ARDMatchLevelCertainty`` values
        indicating how certain the ARD refinement is given typing resolution.
    :ivar molecular_match_levels: Tuple of ``MolecularMatchLevel`` values
        refining ARD-equivalent allele pairs at 1–4-field level
        (``NOT_APPLICABLE`` if not ARD-matched).
    :ivar molecular_match_certainties: Tuple of
        ``MolecularMatchLevelCertainty`` values indicating how certain the
        molecular refinement is given typing resolution.
    :ivar dpb1_tce_status: Optional DPB1 permissive/non-permissive
        classification from the EBI TCE API (only populated for DPB1 loci).
    :ivar is_homozygous_patient: ``True`` if the patient is homozygous at
        this locus at ARD-reduced level, ``False`` if heterozygous, or
        ``None`` if ARD-reduced alleles are not available.
    """
    def __init__(
            self,
            patient: HLAPair,
            donor: HLAPair,
            pairing_score: int,
            allele_match_levels: Tuple[AlleleMatchLevel, AlleleMatchLevel],
            ard_match_levels: Optional[
                Tuple[ARDMatchLevel, ARDMatchLevel]
            ] = None,
            ard_match_level_certainty: Optional[
                Tuple[ARDMatchLevelCertainty, ARDMatchLevelCertainty]
            ] = None,
            molecular_match_levels: Optional[
                Tuple[MolecularMatchLevel, MolecularMatchLevel]
            ] = None,
            molecular_match_level_certainty: Optional[
                Tuple[
                    MolecularMatchLevelCertainty,
                    MolecularMatchLevelCertainty
                ]
            ] = None,
    ) -> None:

        self.patient = patient
        self.donor = donor
        self.pairing_score = pairing_score
        self.allele_match_levels = allele_match_levels
        if ard_match_levels is None:
            self.ard_match_levels = (
                ARDMatchLevel.NOT_APPLICABLE,
                ARDMatchLevel.NOT_APPLICABLE,
            )
        else:
            self.ard_match_levels = ard_match_levels

        if ard_match_level_certainty is None:
            self.ard_match_certainties = (
                ARDMatchLevelCertainty.NOT_APPLICABLE,
                ARDMatchLevelCertainty.NOT_APPLICABLE,
            )
        else:
            self.ard_match_certainties = ard_match_level_certainty
        if molecular_match_levels is None:
            self.molecular_match_levels = (
                MolecularMatchLevel.NOT_APPLICABLE,
                MolecularMatchLevel.NOT_APPLICABLE,
            )
        else:
            self.molecular_match_levels = molecular_match_levels
        if molecular_match_level_certainty is None:
            self.molecular_match_certainties = (
                MolecularMatchLevelCertainty.NOT_APPLICABLE,
                MolecularMatchLevelCertainty.NOT_APPLICABLE,
            )
        else:
            self.molecular_match_certainties = molecular_match_level_certainty

        # optional external matching information
        self.dpb1_tce_status: Optional[DPB1TCEStatus] = None

        # check homozygous patient
        # TODO: homozygosity check currently is capped at ARD which may not
        # be considered *true*
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
        Get locus-level match category for a given resolution.

        :param resolution: Resolution level (``"basic"`` or ``"high"``).
        :return: Match level for the given resolution as a string.
        :raises ValueError: If an unknown resolution level is requested.
        """
        if resolution == "basic":
            return self.loci_match_basic_resolution
        elif resolution == "high":
            return self.loci_match_high_resolution
        else:
            raise ValueError(
                f"Unknown resolution level: {resolution}\n"
                f"Expected 'basic', 'high'."
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

    def _loci_level_match(self, resolution):
        """
        Locus-level match category based on AlleleMatchLevels.
        """
        match_level_1, match_level_2 = self.allele_match_levels

        if (
            match_level_1 is AlleleMatchLevel.NOT_ASSESSABLE
            or match_level_2 is AlleleMatchLevel.NOT_ASSESSABLE
        ):
            return AlleleMatchLevel.NOT_ASSESSABLE.name

        if resolution == 'basic_resolution':
            return self._calculate_loci_match_basic_resolution(
                match_level_1, match_level_2
            )
        elif resolution == 'high_resolution':
            return self._calculate_loci_match_high_resolution(
                match_level_1, match_level_2
            )
        else:
            raise ValueError(
                f"Unknown resolution level: {resolution}\n"
                f"Expected 'basic_resolution', 'high_resolution'."
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
            AlleleMatchLevel.ARD_MATCH
        }
        mismatch_levels = {
            AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH,
            AlleleMatchLevel.ANTIGEN_MISMATCH,
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
            AlleleMatchLevel.ARD_MATCH
        }
        mismatch_levels = {
            AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH,
            AlleleMatchLevel.ANTIGEN_MISMATCH,
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


@dataclass(frozen=True)
class _PairingResult:
    """
    Internal result from allele pairing.

    Intended for research workflows.

    This dataclass stores match levels and certainties of allele pairings.
    Used internally by `_get_correct_allele_pairing`.

    Attributes:
        score (int): Sum of AlleleMatchLevel values for both allele pairs.
            Used as primary criterion for selecting optimal pairing.
        allele_match_levels (Tuple[AlleleMatchLevel, AlleleMatchLevel]):
            ARD-based match level for each paired allele comparison.
        ard_match_levels (Tuple[ARDMatchLevel, ARDMatchLevel]):
            G-group vs P-group refinement for ARD-matched alleles.
            NOT_APPLICABLE if AlleleMatchLevel != ARD_MATCH.
        ard_match_certainties (
            Tuple[ARDMatchLevelCertainty, ARDMatchLevelCertainty]
        ):
            Certainty of ARD match level given typing resolution.
            UNCERTAIN indicates a higher ARDMatchLevel may be possible.
        molecular_match_levels (
            Tuple[MolecularMatchLevel, MolecularMatchLevel]
        ):
            Field-by-field identity refinement for ARD-matched alleles.
            NOT_APPLICABLE if AlleleMatchLevel != ARD_MATCH.
        molecular_match_certainties (
            Tuple[MolecularMatchLevelCertainty, MolecularMatchLevelCertainty]
        ):
            Certainty of molecular match level given typing resolution.
            UNCERTAIN indicates a higher MolecularMatchLevel may be possible.
    """
    score: int
    allele_match_levels: Tuple[AlleleMatchLevel, AlleleMatchLevel]
    ard_match_levels: Tuple[ARDMatchLevel, ARDMatchLevel]
    ard_match_certainties: Tuple[
        ARDMatchLevelCertainty, ARDMatchLevelCertainty
    ]
    molecular_match_levels: Tuple[MolecularMatchLevel, MolecularMatchLevel]
    molecular_match_certainties: Tuple[
        MolecularMatchLevelCertainty, MolecularMatchLevelCertainty
    ]


def _map_expression_decision(
    decision: ExpressionSuffixMatchLevel,
) -> Optional[AlleleMatchLevel]:
    if decision is ExpressionSuffixMatchLevel.IGNORE:
        return None
    mapping = {
        ExpressionSuffixMatchLevel.NOT_ASSESSABLE:
            AlleleMatchLevel.NOT_ASSESSABLE,
        ExpressionSuffixMatchLevel.ALLELE_MISMATCH:
            AlleleMatchLevel.ALLELE_MISMATCH,
        ExpressionSuffixMatchLevel.ANTIGEN_MISMATCH:
            AlleleMatchLevel.ANTIGEN_MISMATCH,
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
    # Any 'Q' present (defaults to NOT_ASSESSABLE)
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

    # (1) LOCUS and LOW-RES comparison

    # first check if loci match (NOTE: DRB3/4/5 hard coded to locus DRB345)
    if hla1.locus != hla2.locus:
        raise InvalidLocusComparisonError(hla1.locus, hla2.locus)

    # for locus DRB345, we still stored the original DRB3/4/5 sub-locus
    if hla1.locus == 'DRB345' and hla1.drb_sub_locus != hla2.drb_sub_locus:
        return AlleleMatchLevel.DRB345_SUBLOCUS_MISMATCH

    if min(
        hla1.has_resolution_level(), hla2.has_resolution_level()
    ) < 1:
        # no allele fields
        return AlleleMatchLevel.NOT_ASSESSABLE

    if min(
        hla1.has_resolution_level(), hla2.has_resolution_level()
    ) < 2:
        # check if allele groups differ
        if hla1.allele_group != hla2.allele_group:
            return AlleleMatchLevel.ANTIGEN_MISMATCH
        # else we cannot determine a match level (missing data)
        else:
            return AlleleMatchLevel.NOT_ASSESSABLE

    # --- from here on we have at least two-field resolution ---
    # (2) TWO-FIELD COMPARISON

    # check for allele group mismatch
    if hla1.allele_group != hla2.allele_group:
        return AlleleMatchLevel.ANTIGEN_MISMATCH

    if (
        hla1.ard_redux_allele_string is None
        or hla2.ard_redux_allele_string is None
    ):
        # NOTE: this should never happen (!)
        raise RuntimeError(
            f"HLA parsing failed for '{hla1.allele_string}' or "
            f"'{hla2.allele_string}'. Please report this issue."
        )

    if hla1.ard_redux_allele != hla2.ard_redux_allele:
        return AlleleMatchLevel.ALLELE_MISMATCH

    # (3) EXPRESSION COMPARISON (suffixes)
    # NOTE: we may need to move expression comparison to be evaluated
    # directly after locus comparison

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

    # (3) ARD MATCH
    return AlleleMatchLevel.ARD_MATCH


def _refine_ard_match_level_by_group_association(
    hla1: HLA,
    hla2: HLA,
    allele_match_level: AlleleMatchLevel
) -> tuple[ARDMatchLevel, ARDMatchLevelCertainty]:
    """
    Compares two ARD-matched HLA alleles and returns an ARDMatchLevel

    Args:
        hla1: First HLA allele object
        hla2: Second HLA allele object
        allele_match_level: AlleleMatchLevel of hla1 and hla2

    Returns:
        Tuple[ARDMatchLevel, ARDMatchLevelCertainty]
            ARDMatchLevel IntEnum value indicating level of ARD matching
            ARDMatchLevelCertainty Enum indicating certainty of ARD match level
    Raises:
        TypeError: If hla1 or hla2 is not an instance of HLA
        InvalidLocusComparisonError: If hla1 and hla2 have incompatible loci

    Only applicable to AlleleMatchLevel == ARD_MATCH
    Otherwise, returns NOT_APPLICABLE for both level and certainty
    """
    # sanity checks
    if not isinstance(hla1, HLA):
        raise TypeError(
            f"hla1 must be an instance of HLA, not {type(hla1).__name__}."
        )
    if not isinstance(hla2, HLA):
        raise TypeError(
            f"hla2 must be an instance of HLA, not {type(hla2).__name__}."
        )
    if hla1.locus != hla2.locus:
        raise InvalidLocusComparisonError(hla1.locus, hla2.locus)

    # additional safeguards against misuse
    is_claimed_ard_match = (allele_match_level is AlleleMatchLevel.ARD_MATCH)
    # HLA class guarantees valid hla two-field allele if redux worked
    has_ard_data = (
        hla1.ard_redux_allele_string is not None and
        hla2.ard_redux_allele_string is not None
    )
    # still redux string must be equal to confirm ARD_MATCH
    is_actual_ard_match = (
        has_ard_data and
        hla1.ard_redux_allele_string == hla2.ard_redux_allele_string
    )

    if is_claimed_ard_match and not has_ard_data:
        raise ARDMatchRefinementError(
            f"ARD_MATCH but ARD reduction data missing. "
            f"hla1.ard_redux_allele_string={hla1.ard_redux_allele_string}, "
            f"hla2.ard_redux_allele_string={hla2.ard_redux_allele_string}"
        )
    if is_claimed_ard_match and not is_actual_ard_match:
        raise ARDMatchRefinementError(
            f"ARD_MATCH but alleles differ at ARD level. "
            f"hla1.ard_redux_allele_string={hla1.ard_redux_allele_string}, "
            f"hla2.ard_redux_allele_string={hla2.ard_redux_allele_string}"
        )
    if not is_claimed_ard_match and is_actual_ard_match:
        raise ARDMatchRefinementError(
            f"{allele_match_level.name} but alleles ARE "
            f"ARD-equivalent. This indicates a bug in the caller. "
            f"hla1.ard_redux_allele_string={hla1.ard_redux_allele_string}, "
            f"hla2.ard_redux_allele_string={hla2.ard_redux_allele_string}"
        )

    # (1) Valid non-ARD_MATCH: return NOT_APPLICABLE
    if allele_match_level is not AlleleMatchLevel.ARD_MATCH:
        return (
            ARDMatchLevel.NOT_APPLICABLE,
            ARDMatchLevelCertainty.NOT_APPLICABLE
        )

    # NOTE: specific group_code ('01P', instead of 'P') encoded in allele field
    # of HLA object (or synonymous_variant field for 'G' group)
    # TODO: imo this needs and update in the HLA parsing logic
    # not a bug per se, just counterintuitive and welcomes errors

    # (2) P-group is first exit if we lack information
    # NOTE: HLA parsing **guarantees** that a given "P" is highest resolution
    if hla1.group_code == "P" or hla2.group_code == "P":
        return (
            ARDMatchLevel.P_GROUP_MATCH,
            # could still be G-group match
            ARDMatchLevelCertainty.UNCERTAIN
        )

    # (3) G-group next
    min_resolution = min(
        hla1.has_resolution_level(), hla2.has_resolution_level()
    )
    # G group is more complex, we have G-group match if:
    # a) hla1.synonymous_variant == hla2.synonymous_variant without G-group
    if min_resolution >= 3:
        if (
            hla1.group_code != "G"
            and hla2.group_code != "G"
            and hla1.synonymous_variant == hla2.synonymous_variant
        ):
            return (
                ARDMatchLevel.G_GROUP_MATCH,
                ARDMatchLevelCertainty.CERTAIN
            )
        # b) hla1.group_code == "G" and hla2.group_code == "G"
        if (
            hla1.group_code == "G"
            and hla2.group_code == "G"
            and hla1.synonymous_variant == hla2.synonymous_variant
        ):
            return (
                ARDMatchLevel.G_GROUP_MATCH,
                ARDMatchLevelCertainty.CERTAIN
            )
        # c) one allele has G-group, the other not, but both are in the same
        # G-group
        if hla1.group_code == "G" or hla2.group_code == "G":
            pyard = get_ard_instance()
            pyard_g1_string = pyard.redux(hla1.allele_string, 'G')
            pyard_g2_string = pyard.redux(hla2.allele_string, 'G')
            if (
                pyard_g1_string == pyard_g2_string
                and pyard_g1_string.endswith('G')
                and pyard_g2_string.endswith('G')
            ):
                return (
                    ARDMatchLevel.G_GROUP_MATCH,
                    ARDMatchLevelCertainty.CERTAIN
                )
        # d) if both are not G-group coded:
        if hla1.group_code != "G" and hla2.group_code != "G":
            pyard = get_ard_instance()
            pyard_g1_string = pyard.redux(hla1.allele_string, 'G')
            pyard_g2_string = pyard.redux(hla2.allele_string, 'G')
            if (
                pyard_g1_string == pyard_g2_string
                and pyard_g1_string.endswith('G')
                and pyard_g2_string.endswith('G')
            ):
                return (
                    ARDMatchLevel.G_GROUP_MATCH,
                    ARDMatchLevelCertainty.CERTAIN
                )
    # (4) quo vadis?
    # due to overlap of P- and G-groups we could actually get more info
    # e.g. A*01:468 and ​A*01:471 are part of A*01:01P and A*01:01:01G
    # however, py-ard's 'G' reduction is currently not robust
    # e.g., print(pyard.redux("A*01:01", 'G')) returns 'A*01:01:01G',
    # but A*01:01:162 (valid allele) is not part of A*01:01:01G

    # NOTE: so until this is resolved for now
    return (
        ARDMatchLevel.P_GROUP_MATCH,
        # could still be G-group match
        ARDMatchLevelCertainty.UNCERTAIN
    )


def _refine_ard_match_level_at_molecular_level(
    hla1: HLA,
    hla2: HLA,
    allele_match_level: AlleleMatchLevel
) -> tuple[MolecularMatchLevel, MolecularMatchLevelCertainty]:
    """
    Compares two ARD-matched HLA alleles and returns a MolecularMatchLevel

    Args:
        hla1: First HLA allele object
        hla2: Second HLA allele object
        allele_match_level: AlleleMatchLevel of hla1 and hla2

    Returns:
        Tuple[MolecularMatchLevel, MolecularMatchLevelCertainty]
            MolecularMatchLevel IntEnum value indicating degree of 1–4 field
            identity
            MolecularMatchLevelCertainty Enum indicating whether a higher level
            is still possible given typing resolution
    Raises:
        TypeError: If hla1 or hla2 is not an instance of HLA
        InvalidLocusComparisonError: If hla1 and hla2 have incompatible loci

    Only applicable to AlleleMatchLevel == ARD_MATCH
    Otherwise, returns NOT_APPLICABLE for both level and certainty
    """
    # sanity checks
    if not isinstance(hla1, HLA):
        raise TypeError(
            f"hla1 must be an instance of HLA, not {type(hla1).__name__}."
        )
    if not isinstance(hla2, HLA):
        raise TypeError(
            f"hla2 must be an instance of HLA, not {type(hla2).__name__}."
        )
    if hla1.locus != hla2.locus:
        raise InvalidLocusComparisonError(hla1.locus, hla2.locus)

    # additional safeguards against misuse
    is_claimed_ard_match = (allele_match_level is AlleleMatchLevel.ARD_MATCH)
    # HLA class guarantees valid hla two-field allele if redux worked
    has_ard_data = (
        hla1.ard_redux_allele_string is not None and
        hla2.ard_redux_allele_string is not None
    )
    # still redux string must be equal to confirm ARD_MATCH
    is_actual_ard_match = (
        has_ard_data and
        hla1.ard_redux_allele_string == hla2.ard_redux_allele_string
    )

    if is_claimed_ard_match and not has_ard_data:
        raise ARDMatchRefinementError(
            f"ARD_MATCH but ARD reduction data missing. "
            f"hla1.ard_redux_allele_string={hla1.ard_redux_allele_string}, "
            f"hla2.ard_redux_allele_string={hla2.ard_redux_allele_string}"
        )
    if is_claimed_ard_match and not is_actual_ard_match:
        raise ARDMatchRefinementError(
            f"ARD_MATCH but alleles differ at ARD level. "
            f"hla1.ard_redux_allele_string={hla1.ard_redux_allele_string}, "
            f"hla2.ard_redux_allele_string={hla2.ard_redux_allele_string}"
        )
    if not is_claimed_ard_match and is_actual_ard_match:
        raise ARDMatchRefinementError(
            f"{allele_match_level.name} but alleles ARE "
            f"ARD-equivalent. This indicates a bug in the caller. "
            f"hla1.ard_redux_allele_string={hla1.ard_redux_allele_string}, "
            f"hla2.ard_redux_allele_string={hla2.ard_redux_allele_string}"
        )

    # (1) Valid non-ARD_MATCH: return NOT_APPLICABLE
    # e.g. A*01:01 vs A*02:01
    if allele_match_level is not AlleleMatchLevel.ARD_MATCH:
        return (
            MolecularMatchLevel.NOT_APPLICABLE,
            MolecularMatchLevelCertainty.NOT_APPLICABLE
        )

    # (2) P-group code: molecular not applicable
    # e.g. A*01:01P vs A*01:01:01:01
    if hla1.group_code == "P" or hla2.group_code == "P":
        return (
            MolecularMatchLevel.NOT_ASSESSABLE,
            # Could be protein/coding/exact
            MolecularMatchLevelCertainty.UNCERTAIN
        )

    # From here: ARD_MATCH, no P-group - let's try to refine ARD_MATCH
    # we need the resolution multiple times
    min_resolution = min(
        hla1.has_resolution_level(), hla2.has_resolution_level()
    )

    # (3) cases with res == 2
    # TODO: double check if G group interferes with protein
    # we either have full protein match:
    # e.g. A*01:01 vs A*01:01
    if (min_resolution == 2) and (hla1.allele == hla2.allele):
        return (
            MolecularMatchLevel.FULL_PROTEIN_MATCH,
            # Could be coding/exact
            MolecularMatchLevelCertainty.UNCERTAIN
        )
    # or a mismatch:
    # e.g. A*01:01 vs A*01:15 (same P-group)
    elif (min_resolution == 2) and (hla1.allele != hla2.allele):
        return (
            MolecularMatchLevel.ARD_MATCH_ONLY,
            MolecularMatchLevelCertainty.CERTAIN
        )

    # (4) cases with res >= 3
    # NOTE: we will look for a more elegant way to structure this logic
    if min_resolution >= 3:
        # If second fields differ, we can only say ARD_MATCH_ONLY:
        # they are ARD-equivalent but not protein/coding/exact identical
        # e.g. A*02:01:01 vs A*02:09:01 (both A*02:01P)
        if hla1.allele != hla2.allele:
            return (
                MolecularMatchLevel.ARD_MATCH_ONLY,
                MolecularMatchLevelCertainty.CERTAIN
            )
        # NOTE: this should also handle all G-groups

        # From here: first two fields identical -> at least FULL_PROTEIN_MATCH

        # (4a) min_resolution == 3: we know the 3rd field (synonymous variant)
        if min_resolution == 3:
            if (
                hla1.synonymous_variant == hla2.synonymous_variant
                and hla1.group_code != "G"
                and hla2.group_code != "G"
            ):
                # 1–3 fields identical and 4th is unknown or untyped
                return (
                    MolecularMatchLevel.CODING_SEQUENCE_MATCH,
                    # Could still be EXACT_ALLELE_MATCH if 4th field also equal
                    MolecularMatchLevelCertainty.UNCERTAIN
                )
            if (
                hla1.synonymous_variant == hla2.synonymous_variant
                and hla1.group_code == "G"
                and hla2.group_code == "G"
            ):
                # 1–3 fields identical but only G-group known
                return (
                    MolecularMatchLevel.FULL_PROTEIN_MATCH,
                    # Could still be coding/exact if 3rd-4th field also equal
                    MolecularMatchLevelCertainty.UNCERTAIN
                )
            elif (
                hla1.synonymous_variant == hla2.synonymous_variant
                # if we have a single 'G' group, we cannot be sure about coding
                and (hla1.group_code == "G" or hla2.group_code == "G")
            ):
                # Same protein but different coding sequence
                return (
                    MolecularMatchLevel.FULL_PROTEIN_MATCH,
                    # Cannot be CODING_SEQUENCE_MATCH or EXACT_ALLELE_MATCH
                    MolecularMatchLevelCertainty.UNCERTAIN
                )
            else:
                # 1–2 fields identical, 3rd differs
                return (
                    MolecularMatchLevel.FULL_PROTEIN_MATCH,
                    MolecularMatchLevelCertainty.CERTAIN
                )

        # (4b) min_resolution == 4: both alleles have 4-field resolution
        # First two fields already identical -> inspect 3rd and 4th:

        if hla1.synonymous_variant != hla2.synonymous_variant:
            # Different third field -> different coding sequence
            return (
                MolecularMatchLevel.FULL_PROTEIN_MATCH,
                MolecularMatchLevelCertainty.CERTAIN
            )

        # Third field identical -> check non-coding (4th) field
        if hla1.non_coding_variant == hla2.non_coding_variant:
            # All 1–4 fields identical
            return (
                MolecularMatchLevel.EXACT_ALLELE_MATCH,
                MolecularMatchLevelCertainty.CERTAIN
            )
        else:
            # 1–3 fields identical, 4th differs
            return (
                MolecularMatchLevel.CODING_SEQUENCE_MATCH,
                MolecularMatchLevelCertainty.CERTAIN
            )

    # if we are here, imho something is odd
    raise ARDMatchRefinementError(
        "_refine_ard_match_level_at_molecular_level was unable to process "
        f"hla1='{hla1.allele_string}', hla2='{hla2.allele_string}'."
    )


def _get_correct_allele_pairing(
        patient_alleles: HLAPair, donor_alleles: HLAPair
) -> _PairingResult:
    """
    Determines the correct pairing of patient and donor HLA alleles by
    evaluating all possible combinations.

    Intended for research workflows.

    :param patient_alleles: ``HLAPair`` containing two patient HLA alleles.
    :param donor_alleles: ``HLAPair`` containing two donor HLA alleles.
    :return: A ``_PairingResult`` instance containing allele-, ARD- and
        molecular-level match classifications and certainties for the
        optimal pairing.

    Notes:
        - Considers two possible pairings:
            1. (patient_hla1, donor_hla1) and (patient_hla2, donor_hla2)
            2. (patient_hla1, donor_hla2) and (patient_hla2, donor_hla1)
        - Selection uses three-level lexicographic scoring:
            1. Primary: AlleleMatchLevel sum (match vs mismatch)
            2. Secondary: MolecularMatchLevel sum (field identity)
            3. Tertiary: ARDMatchLevel sum (G-group vs P-group in ARD)
        - If all scores are equal, the first pairing is returned

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

    # Lexicographic comparison: (allele, molecular, ard)
    best_score: Tuple[float, float, float] = (
        float('-inf'), float('-inf'), float('-inf')
    )
    best_result: Optional[_PairingResult] = None

    for pairing in pairings:
        patient_hla1, donor_hla1, patient_hla2, donor_hla2 = pairing

        # (1) Primary: AlleleMatchLevel
        allele_match1 = allele_match(patient_hla1, donor_hla1)
        allele_match2 = allele_match(patient_hla2, donor_hla2)
        allele_score = int(allele_match1) + int(allele_match2)

        # (2) Refinements if ARD_MATCH
        if allele_match1 is AlleleMatchLevel.ARD_MATCH:
            ard_match1, ard_certainty1 = \
                _refine_ard_match_level_by_group_association(
                    patient_hla1, donor_hla1, allele_match1
                )
            molecular_match1, molecular_certainty1 = \
                _refine_ard_match_level_at_molecular_level(
                    patient_hla1, donor_hla1, allele_match1
                )
        else:
            ard_match1 = ARDMatchLevel.NOT_APPLICABLE
            ard_certainty1 = ARDMatchLevelCertainty.NOT_APPLICABLE
            molecular_match1 = MolecularMatchLevel.NOT_APPLICABLE
            molecular_certainty1 = MolecularMatchLevelCertainty.NOT_APPLICABLE

        # Same allele_match2
        if allele_match2 is AlleleMatchLevel.ARD_MATCH:
            ard_match2, ard_certainty2 = \
                _refine_ard_match_level_by_group_association(
                    patient_hla2, donor_hla2, allele_match2
                )
            molecular_match2, molecular_certainty2 = \
                _refine_ard_match_level_at_molecular_level(
                    patient_hla2, donor_hla2, allele_match2
                )
        else:
            ard_match2 = ARDMatchLevel.NOT_APPLICABLE
            ard_certainty2 = ARDMatchLevelCertainty.NOT_APPLICABLE
            molecular_match2 = MolecularMatchLevel.NOT_APPLICABLE
            molecular_certainty2 = MolecularMatchLevelCertainty.NOT_APPLICABLE

        # (3) Tie-breaker scores
        molecular_score = int(molecular_match1) + int(molecular_match2)
        ard_score = int(ard_match1) + int(ard_match2)

        # (4) Lexicographic comparison: (allele, molecular, ard)
        current_score = (allele_score, molecular_score, ard_score)

        if current_score > best_score:
            best_score = current_score
            best_result = _PairingResult(
                score=allele_score,
                allele_match_levels=(allele_match1, allele_match2),
                ard_match_levels=(ard_match1, ard_match2),
                ard_match_certainties=(ard_certainty1, ard_certainty2),
                molecular_match_levels=(molecular_match1, molecular_match2),
                molecular_match_certainties=(
                    molecular_certainty1, molecular_certainty2
                ),
            )

    # Cannot be None: we always have exactly 2 pairings
    return best_result  # type: ignore[return-value]


def allele_pair_match(patient: HLAPair, donor: HLAPair) -> MatchResult:
    """
    Compute research match/mismatch levels for two HLA allele pairs, one in
    the 'patient' role and one in the 'donor' role.

    Intended for research workflows.

    :param patient: Patient ``HLAPair`` containing two HLA alleles.
    :param donor: Donor ``HLAPair`` containing two HLA alleles.
    :return: ``MatchResult`` object storing allele-level match categories
        and all ARD and molecular refinements for the optimal pairing.

    Notes:
        - The function assumes that both patient and donor have exactly two

          HLA alleles
        - Uses `_get_correct_allele_pairing` to evaluate all possible

          allele pairings and selects the one with the highest score
    """
    result = _get_correct_allele_pairing(patient, donor)

    return MatchResult(
        patient=patient,
        donor=donor,
        pairing_score=result.score,
        allele_match_levels=result.allele_match_levels,
        ard_match_levels=result.ard_match_levels,
        ard_match_level_certainty=result.ard_match_certainties,
        molecular_match_levels=result.molecular_match_levels,
        molecular_match_level_certainty=result.molecular_match_certainties,
    )


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
                "matching will be reported as NOT_ASSESSABLE."
            )
            if locus == "DRB345":
                donor_pair = HLAPair(HLA("DRBX*NA"), HLA("DRBX*NA"))
            else:
                donor_pair = HLAPair(HLA(f"{locus}*NA"), HLA(f"{locus}*NA"))

        # compute match (missing donor pair will propagate NOT_ASSESSABLE)
        match_result = allele_pair_match(patient_pair, donor_pair)

        # additional diagnostics
        if all(level == AlleleMatchLevel.NOT_ASSESSABLE
               for level in match_result.allele_match_levels):
            logger.warning(
                f"Typing resolution insufficient for locus {locus} "
                f"(patient {patient_pair} / donor {donor_pair})."
            )
        results.append(match_result)

    return results
