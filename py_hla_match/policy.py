from __future__ import annotations


from dataclasses import dataclass
from enum import IntEnum, Enum
from typing import FrozenSet


class AlleleMatchLevel(IntEnum):
    """
    Following hla nomenclature:

    NOT_ASSESSABLE (0):
        Typing resolution insufficient
    LOCUS_MISMATCH (-3):
        Mismatch at HLA locus DRB3/-4/-5
    ANTIGEN_MISMATCH (-2):
        Mismatch at the group code encoding antigen
    ALLELE_MISMATCH (-1):
        Mismatch at the allele level encoding a specific allele
    ARD_MATCH (1):
        ARD level match

    cf.https://hla.alleles.org/nomenclature/naming.html
    """
    NOT_ASSESSABLE = 0
    # Mismatch levels
    LOCUS_MISMATCH = -3
    ANTIGEN_MISMATCH = -2
    ALLELE_MISMATCH = -1
    # Separating mismatches from matches based on ARD
    ARD_MATCH = 1


class ARDMatchLevel(IntEnum):
    """
    ARD refinement iff AlleleMatchLevel == ARD_MATCH.

    G_GROUP_MATCH (2):
        ARD exons identical at nucleotide level (G group)
    P_GROUP_MATCH (1):
        ARD exons identical at amino-acid level (P group)
    NOT_APPLICABLE:
        AlleleMatchLevel != ARD_MATCH
    """
    G_GROUP_MATCH = 2
    P_GROUP_MATCH = 1
    NOT_APPLICABLE = 0


class ARDMatchLevelCertainty(Enum):
    """
    Certainty of ARD-match level.

    NOT_APPLICABLE:
        AlleleMatchLevel != ARD_MATCH
    UNCERTAIN:
        insufficient typing resolution to be sure about ARDMatchLevel
        i.e., a *higher* ARDMatchLevel **is possible**
    CERTAIN:
        sufficient typing resolution to be sure about ARDMatchLevel
        i.e., a *higher* ARDMatchLevel **is not possible**
    """
    NOT_APPLICABLE = "not applicable for AlleleMatchLevel != ARD_MATCH"
    UNCERTAIN = "uncertain about ARDMatchLevel due to insufficient resolution"
    CERTAIN = "certain about ARDMatchLevel due to sufficient resolution"


class MolecularMatchLevel(IntEnum):
    """
    Molecular (sequence-level) refinement iff AlleleMatchLevel == ARD_MATCH.

    EXACT_ALLELE_MATCH (5):
        1-4 fields identical
    CODING_SEQUENCE_MATCH (4):
        1-3 fields identical, different or untyped 4th field
    FULL_PROTEIN_MATCH (3):
        1-2 fields identical, different or untyped 3rd field
    ARD_MATCH_ONLY (2):
        2-field difference but ARD equivalent
    NOT_ASSESSABLE (1):
        Typing resolution insufficient to assess molecular match level
    NOT_APPLICABLE (0):
        AlleleMatchLevel != ARD_MATCH
    """
    EXACT_ALLELE_MATCH = 5
    CODING_SEQUENCE_MATCH = 4
    FULL_PROTEIN_MATCH = 3
    ARD_MATCH_ONLY = 2
    NOT_ASSESSABLE = 1
    NOT_APPLICABLE = 0


class MolecularMatchLevelCertainty(Enum):
    """
    Certainty of molecular level.

    NOT_APPLICABLE:
        AlleleMatchLevel != ARD_MATCH
    UNCERTAIN:
        insufficient typing resolution to be sure about MolecularMatchLevel
        i.e., a *higher* MolecularMatchLevel **is possible**
    CERTAIN:
        sufficient typing resolution to be sure about MolecularMatchLevel
        i.e., a *higher* MolecularMatchLevel **is not possible**
    """
    NOT_APPLICABLE = "not applicable for AlleleMatchLevel != ARD_MATCH"
    UNCERTAIN = "uncertain about MolecularMatchLevel due to insufficient " \
        "resolution"
    CERTAIN = "certain about MolecularMatchLevel due to sufficient resolution"


class ExpressionSuffixMatchLevel(Enum):
    """
    Policy decisions about how to treat expression suffixes in matching.
    Similar to cf AlleleMatchLevel in hla.py, but for expression suffixes.
    """
    NOT_ASSESSABLE = "not_assessable"
    ALLELE_MISMATCH = "allele_mismatch"
    ANTIGEN_MISMATCH = "antigen_mismatch"
    ARD_MATCH = "ard_match"
    IGNORE = "ignore"


@dataclass(frozen=True)
class ExpressionSuffixPolicy:
    # N null; L low; S secreted; C cytoplasmic; A aberrant; Q questionable
    risk_suffixes: FrozenSet[str] = frozenset({"N", "L", "S", "C", "A"})
    ambiguous_suffixes: FrozenSet[str] = frozenset({"Q"})

    # Defaults
    equal_risk: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.NOT_ASSESSABLE
    risk_vs_none: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.ALLELE_MISMATCH
    risk_vs_different_risk: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.ALLELE_MISMATCH
    q_present: ExpressionSuffixMatchLevel = \
        ExpressionSuffixMatchLevel.NOT_ASSESSABLE
