from __future__ import annotations


from dataclasses import dataclass
from enum import IntEnum, Enum
from typing import FrozenSet


class AlleleMatchLevel(IntEnum):
    """
    Following hla nomenclature:
    LOCUS_MISMATCH: Mismatch at a particular HLA locus
    ALLELE_GROUP_MISMATCH: Mismatch at the group code
    ALLELE_MISMATCH: Mismatch at the allele level
    ARD_MATCH: ARD level match
    SYNONYMOUS_VARIANT_MATCH: Synonymous variant match
    NON_CODING_VARIANT_MATCH: Non-coding variant match
    cf.https://hla.alleles.org/nomenclature/naming.html
    """
    # NOTE: no clean seperation for now between match codes and not applicable
    NOT_APPLICABLE = 0

    # clean AlleleMatchLevel
    LOCUS_MISMATCH = -3
    ALLELE_GROUP_MISMATCH = -2
    ALLELE_MISMATCH = -1
    # NOTE: seperating mismatches from matches based on ARD
    ARD_MATCH = 1
    SYNONYMOUS_VARIANT_MATCH = 2
    NON_CODING_VARIANT_MATCH = 3


class ExpressionSuffixMatchLevel(Enum):
    """
    Policy decisions about how to treat expression suffixes in matching.
    Similar to cf AlleleMatchLevel in hla.py, but for expression suffixes.
    """
    NOT_APPLICABLE = "not_applicable"
    ALLELE_MISMATCH = "allele_mismatch"
    ALLELE_GROUP_MISMATCH = "allele_group_mismatch"
    ARD_MATCH = "ard_match"
    IGNORE = "ignore"


@dataclass(frozen=True)
class ExpressionSuffixPolicy:
    # N null; L low; S secreted; C cytoplasmic; A aberrant; Q questionable
    risk_suffixes: FrozenSet[str] = frozenset({"N", "L", "S", "C", "A"})
    ambiguous_suffixes: FrozenSet[str] = frozenset({"Q"})

    # Defaults
    equal_risk: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.NOT_APPLICABLE
    risk_vs_none: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.ALLELE_MISMATCH
    risk_vs_different_risk: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.ALLELE_MISMATCH
    q_present: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.NOT_APPLICABLE
