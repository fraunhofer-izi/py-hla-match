# hla.py
import logging
import re
from .singleton import get_ard_instance

from py_hla_match.exceptions import MalformedHLAStringError

logger = logging.getLogger(__name__)

# regex pattern for HLA allele string
NOMENCLATURE_PATTERN = re.compile(
    r"""
    ^(?:HLA-)?
    (?P<locus>[A-Z0-9]+)
    \*
    (?P<allele_fields>\d{2,}(?::\d{2,}){0,3})
    (?P<suffix>[NLSCAQ])?
    (?P<group_code>[GP])?$
    """,
    re.VERBOSE
)

LOCUS_PATTERN = re.compile(
    r"""
    ^(?:HLA-)?
    (?P<locus>[A-Z0-9]+)
    \*
    (?P<remainder>.*)?$
    """,
    re.VERBOSE
)

# TODO: discuss additional regex patterns

# regex pattern for ARD redux allele string
REDUX_PATTERN = re.compile(
    r"""
    ^(?:HLA-)?                 # Optional 'HLA-' prefix
    (?P<locus>[A-Z0-9]+)       # Locus
    \*                         # Asterisk separator
    (?P<allele_fields>\d{2,}(?::\d{2,}){0,3})  # Allele fields
    """,
    re.VERBOSE
)

VALID_HLA_LOCI = frozenset({
    'A', 'B', 'C', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'N', 'P', 'R', 'S', 'T',
    'U', 'V', 'W', 'X', 'Y', 'Z', 'DRA', 'DRB1', 'DRB2', 'DRB3', 'DRB4',
    'DRB5', 'DRB6', 'DRB7', 'DRB8', 'DRB9', 'DRBX',
    'DQA1', 'DQB1', 'DQA2', 'DQB2', 'DQB3',
    'DOA', 'DOB', 'DMA', 'DMB',
    'DPA1', 'DPB1', 'DPA2', 'DPB2', 'DPA3',
    'HFE', 'TAP1', 'TAP2', 'PSMB9', 'PSMB8',
    'MICA', 'MICB', 'MICC', 'MICD', 'MICE'
})


class HLA:
    def __init__(self, allele_string: str) -> None:
        """
        Initializes an HLA object by parsing an HLA allele string.

        :param allele_string: HLA allele string
        """
        self.allele_string = allele_string
        self.locus = None
        self.allele_group = None
        self.allele = None
        self.synonymous_variant = None
        self.non_coding_variant = None
        self.suffix = None
        self.group_code = None
        self.ard_redux_allele_string = None
        self.ard_redux_allele_group = None
        self.ard_redux_allele = None

        # NOTE: currently the parser allows invalid allele strings at
        # allele_group
        self._parse_allele()
        if self.allele:
            self._ard_redux()
        else:
            logger.warning(
                f"HLA string '{self.allele_string}' at locus "
                f"'{self.locus}' is not a specific allele."
            )
            if self.allele_group is not None:
                logger.warning(
                    f"WARNING: Validity of '{self.allele_group}' is not "
                    "checked."
                )

    def _parse_allele(self) -> None:
        """
        Parses HLA allele string and populate HLA attributes.
        """
        # validate the allele string
        match = self._validate_nomenclature()

        # extract locuse
        self.locus = match.group('locus')

        try:
            # try NOMENCLATURE_PATTERN first
            allele_fields = match.group('allele_fields')
            if allele_fields is not None:
                self.suffix = match.group('suffix')
                self.group_code = match.group('group_code')

                # extract details from allele fields
                field_contents = allele_fields.split(':')
                if len(field_contents) > 0:
                    self.allele_group = field_contents[0]
                if len(field_contents) > 1:
                    self.allele = field_contents[1]
                if len(field_contents) > 2:
                    self.synonymous_variant = field_contents[2]
                if len(field_contents) > 3:
                    self.non_coding_variant = field_contents[3]
        except (AttributeError, IndexError):  # thrown by regex match
            # LOCUS_PATTERN, check for remainder
            logger.warning(
                f"HLA string '{self.allele_string}' did match HLA Nomenclature"
            )
            remainder = match.group('remainder')
            if remainder:
                logger.warning(
                    f"HLA string '{self.allele_string}' at locus"
                    f" '{self.locus}' has unparseable content: '{remainder}'"
                )
            else:
                logger.warning(
                    f"HLA string '{self.allele_string}' at locus"
                    f" '{self.locus}' is not a valid HLA allele."
                )

    def _is_valid_locus(self, locus: str) -> bool:
        """Locus validation using known loci."""
        return locus in VALID_HLA_LOCI

    def _validate_nomenclature(self) -> re.Match:
        """
        Validate HLA allele string with nomenclature.

        Tries to extract locus information if HLA allele string is not
        complete.

        :raises MalformedHLAStringError: If the allele string is malformed.
        """
        # check if we got at least a one-field (valid) allele string
        allele_match = NOMENCLATURE_PATTERN.match(self.allele_string)
        if allele_match is not None:
            locus = allele_match.group('locus')
            if not self._is_valid_locus(locus):
                # invalid locus
                raise MalformedHLAStringError(
                    f"Invalid HLA allele string: {self.allele_string} with "
                    f"unknown locus '{locus}'."
                )
            return allele_match

        # try to exctract valid locus
        locus_match = LOCUS_PATTERN.match(self.allele_string)
        if locus_match is not None:
            locus = locus_match.group('locus')
            if not self._is_valid_locus(locus):
                # invalid locus
                raise MalformedHLAStringError(
                    f"Invalid HLA allele string: {self.allele_string} with "
                    f"unknown locus: {locus}"
                )
            return locus_match

        # if we reach here, we have a malformed HLA string
        raise MalformedHLAStringError(
            f"Invalid HLA allele string: {self.allele_string} with "
            f"unknown locus: {self.locus}"
        )

        # TODO: Return additional details on specific error occurence
        # TODO: potentially validate the allele string is known allele with
        # ref. to https://hla.alleles.org/alleles/index.html)
        # this is currently handled by py-ard subsequently

    def _ard_redux(self):
        """
        Reduce Allele to 2 field ARD level using the wonderful py-ard package
        For further information: https://github.com/nmdp-bioinformatics/py-ard
        """
        ard = get_ard_instance()
        redux_type = 'lgx'
        self.ard_redux_allele_string = ard.redux(
            self.allele_string, redux_type
        )
        match = REDUX_PATTERN.match(self.ard_redux_allele_string)
        if match:
            allele_fields = match.group('allele_fields')
            field_contents = allele_fields.split(':')
            if len(field_contents) > 0:
                self.ard_redux_allele_group = field_contents[0]
            if len(field_contents) > 1:
                self.ard_redux_allele = field_contents[1]

    def has_two_field_resolution(self) -> bool:
        """Check if HLA has at least two-field resolution."""
        return self.allele is not None

    def has_one_field_resolution(self) -> bool:
        """Check if HLA has at least one-field resolution."""
        return self.allele_group is not None

    def is_specific_allele(self) -> bool:
        """Check if HLA is a specific HLA allele."""
        return self.has_two_field_resolution()

    def __eq__(self, other):
        if not isinstance(other, HLA):
            return NotImplemented
        return (
            self.locus == other.locus and
            self.allele_group == other.allele_group and
            self.allele == other.allele and
            self.synonymous_variant == other.synonymous_variant and
            self.non_coding_variant == other.non_coding_variant and
            self.suffix == other.suffix and
            self.group_code == other.group_code
            # TODO: or (self.reduced = other.reduced )
        )

    def __hash__(self):
        return hash((
            self.locus,
            self.allele_group,
            self.allele,
            self.synonymous_variant,
            self.non_coding_variant,
            self.suffix,
            self.group_code
        ))

    def __repr__(self):
        return (
            f"HLA(allele_string={repr(self.allele_string)}, "
            f"locus={repr(self.locus)}, "
            f"allele_group={repr(self.allele_group)}, "
            f"allele={repr(self.allele)}, "
            f"synonymous_variant={repr(self.synonymous_variant)}, "
            f"non_coding_variant={repr(self.non_coding_variant)}, "
            f"suffix={repr(self.suffix)}, "
            f"group_code={repr(self.group_code)}, "
            f"ard_redux_allele_string={repr(self.ard_redux_allele_string)}, "
            f"ard_redux_allele_group={repr(self.ard_redux_allele_group)}, "
            f"ard_redux_allele={repr(self.ard_redux_allele)})"
        )
