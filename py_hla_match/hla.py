import logging
import re
from .singleton import get_ard_instance

from py_hla_match.exceptions import MalformedHLAStringError

logger = logging.getLogger(__name__)

# regex pattern for HLA allele string
NOMENCLATURE_PATTERN = re.compile(
    r"""
    ^ (?:HLA-)?
    (?P<locus>[A-Z0-9]+) # locus always required
    \* # asterisk always required
    ( # three alternatives
        # 1: well-formed allele
        (?P<allele_fields>\d{2,4}(?::\d{2,4}){0,3})
        (?![NLSCAQ][GP]|[GP][NLSCAQ])  # no suffix AND group code
        (?P<suffix>[NLSCAQ])?
        (?P<group_code>[GP])?
        $ # must match to end of string
    |
        # 2: known nan
        (?P<nan>NA|NE|NEW|UNKNOWN|ND|NULL)
        $
    |
        # 3: anything else - trigger MalformedHLAStringError
        (?P<remainder>.*) # anything else
        $ # must match to end of string
    )
    """,
    re.VERBOSE
)

# regex pattern for ARD redux allele string
REDUX_PATTERN = re.compile(
    r"""
    ^ (?:HLA-)?
    (?P<locus>[A-Z0-9]+)
    \*
    (
        (?P<allele_fields>\d{2,4}(?::\d{2,4}){0,3})
        (?P<suffix>[NLSCAQ])?
        (?P<group_code>[GP])?
        $
    |
        (?P<remainder>.*)
        $
    )
    """,
    re.VERBOSE
)

VALID_HLA_LOCI = frozenset({
    'A', 'B', 'C',  # NOTE: standard HLA Class I loci
    'DRB1',  # NOTE: standard HLA Class II locus
    'DQB1',  # NOTE: standard HLA Class II locus
    'DPB1',  # NOTE: not standard but often discussed
    'DQA1',  # NOTE: interesting candidate
    'DPA1',  # NOTE: interesting candidate
    'DRB3', 'DRB4', 'DRB5',  # NOTE: interesting candidate(s)
    'DRB345',  # NOTE: NOT A LOCUS, but used here as generic DRB3/4/5
    'DRBX',  # NOTE: NOT A LOCUS, but often used to indicate missing DRB3/4/5
    # However, DRB3/4/5 may need special hadling

    'MICA', 'MICB',  # NOTE: interesting candidates
    'DMA', 'DMB', 'DOA', 'DOB',  # NOTE: interesting candidates
    # the exhaustive list of all known HLA loci
    # should be considered non-standard or experimental

    'E', 'F', 'G',  # limited data
    'DQA2', 'DQB2',  # limited data

    # currently not considered for matching ---

    # 'MICC', 'MICD', 'MICE'  # non functional related to MICA/MICB
    # 'TAP1', 'TAP2',  # not direct cell-surface targets
    # 'PSMB9', 'PSMB8',  # not direct cell-surface targets
    # 'HFE',  # hemochromatosis locus, not relevant to HLA matching
    # 'DRA',  # largely monomorphic/minimally polymorphic DR alpha chain
    # likely pseudogenes ---
    # 'DPA2', 'DPB2', 'DPA3', 'DQB3',
    # --- likely pseudogenes
    # pseudogenes ---
    # 'DRB2', 'DRB6', 'DRB7', 'DRB8', 'DRB9',
    # 'H', 'J', 'K', 'L', 'N', 'P', 'R', 'S', 'T',
    # 'U', 'V', 'W', 'X', 'Y', 'Z',
    # --- pseudogenes

    # --- currently not considered for matching
})

# we now support processing of DRB345
DRB345_SUB_LOCI = {
    'DRB3', 'DRB4', 'DRB5',
    'DRBX'  # generic locus indicating "missing"
}


class HLA:
    __slots__ = (
        # original allele string
        'allele_string',
        # parsed locus fields
        'locus', 'drb_sub_locus',
        # parsed allele fields
        'allele_group', 'allele', 'synonymous_variant', 'non_coding_variant',
        # parsed suffix and group code
        'suffix', 'group_code',
        # ARD reduction string
        'ard_redux_allele_string',
        # ARD reduction fields
        'ard_redux_allele_group', 'ard_redux_allele',
        # locked state
        '_locked',
    )

    def __setattr__(self, name, value):
        """
        Block any attribute mutation once _locked is True.
        __slots__ prevents creation of *new* attributes;
        this method prevents changing existing ones.
        """
        if getattr(self, '_locked', False):
            raise AttributeError(
                f"{self.__class__.__name__} instances are immutable "
                f"after initialisation (attempted to set '{name}')."
            )
        super().__setattr__(name, value)

    def __delattr__(self, name):
        raise AttributeError(
            f"{self.__class__.__name__} instances are immutable"
        )

    def __init__(self, allele_string: str) -> None:
        super().__setattr__('allele_string', allele_string)
        # set None as default for missing slot members
        for attr in self.__slots__:
            if attr not in ('allele_string', '_locked'):
                super().__setattr__(attr, None)

        self._parse_allele_string()

        # if well-formed high-res allele use py-ard reduction
        if self.allele:
            self._ard_redux()
        else:
            logger.warning(
                f"HLA string '{self.allele_string}' at locus "
                f"'{self.locus}' is not a specific allele."
            )
            if self.allele_group is not None:
                logger.warning(
                    f"WARNING: Validity of '{self.allele_group}' not checked."
                )

        super().__setattr__('_locked', True)

    def _parse_allele_string(self) -> None:
        """
        Parses HLA allele string and populate HLA attributes.
        """
        # validate the allele string
        match = self._validate_nomenclature()

        # extract locuse
        self.locus = match.group('locus')
        # handle DRB3/4/5 region
        if self.locus in DRB345_SUB_LOCI:
            # if locus is DRB3/4/5, set the sub-locus
            self.drb_sub_locus = self.locus
            self.locus = 'DRB345'

        allele_fields = match.group('allele_fields')
        nan_field = match.group('nan')
        remainder = match.group('remainder')

        if allele_fields:
            self.suffix = match.group('suffix')
            self.group_code = match.group('group_code')
            if self.group_code == 'G' and allele_fields.count(':') < 2:
                raise MalformedHLAStringError(
                    f"'{self.allele_string}' – 'G' group needs ≥3 fields."
                )
            if self.group_code == 'P' and allele_fields.count(':') < 1:
                raise MalformedHLAStringError(
                    f"'{self.allele_string}' – 'P' group needs ≥2 fields."
                )
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
        elif nan_field:
            logger.info(
                f"HLA string '{self.allele_string}' at locus "
                f"'{self.locus}' is undefined: '{nan_field}')."
            )
        elif remainder:
            raise MalformedHLAStringError(
                f"Invalid HLA allele string: '{self.allele_string}' "
                f"contains unparsable content: '{remainder}'"
            )
        else:
            raise MalformedHLAStringError(
                f"HLA string '{self.allele_string}' at locus"
                f" '{self.locus}' is empty."
            )

    def _validate_nomenclature(self) -> re.Match:
        """
        Validate HLA allele string with nomenclature.

        Tries to extract locus information if HLA allele string is not
        complete.

        :raises MalformedHLAStringError: If the allele string is malformed.
        """
        # may wanna check with:
        # https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt
        # check if we got at least a one-field (valid) allele string

        # NOTE: we could soften this in the future
        if any(ch.islower() for ch in self.allele_string if ch.isalpha()):
            raise MalformedHLAStringError(
                f"Lower-case letters found in '{self.allele_string}'. "
                "Allele strings must be upper-case adhering to HLA "
                "nomenclature."
            )

        # NOTE: we might soften this in the future
        if allele_string := self.allele_string:
            if (
                allele_string != allele_string.strip() 
                or any(ch.isspace() for ch in allele_string)
            ):
                raise MalformedHLAStringError(
                    "Found whitespace or invisible characters in allele "
                    f"string '{allele_string}'."
                )

        match = NOMENCLATURE_PATTERN.match(self.allele_string)
        if not match:
            raise MalformedHLAStringError(
                f"Invalid HLA allele string: {self.allele_string}"
                " String must contain valid LOCUS followed by '*'."
            )

        # validate locus
        locus = match.group('locus')
        if not self._is_valid_locus(locus):
            raise MalformedHLAStringError(
                f"Invalid HLA allele string: '{self.allele_string}' with "
                f"unknown locus '{locus}'."
            )

        return match

    def _is_valid_locus(self, locus: str) -> bool:
        """Locus validation using known loci."""
        return locus in VALID_HLA_LOCI

    def _ard_redux(self):
        """
        Reduce Allele to 2 field ARD level using the wonderful py-ard package
        For further information: https://github.com/nmdp-bioinformatics/py-ard
        """
        ard = get_ard_instance()
        redux_type = 'lgx'
        self.ard_redux_allele_string = ard.redux(
            self.allele_string, redux_type
        ).strip()
        match = REDUX_PATTERN.match(self.ard_redux_allele_string)
        if not match:
            raise MalformedHLAStringError(
                "py-ard reports unexpected string not matching regex "
                f"'{self.ard_redux_allele_string}'."
            )
        if match:
            allele_fields = match.group('allele_fields')
            field_contents = allele_fields.split(':')
            if len(field_contents) > 0:
                self.ard_redux_allele_group = field_contents[0]
            if len(field_contents) > 1:
                self.ard_redux_allele = field_contents[1]

    def has_resolution_level(self) -> int:
        """
        Returns the resolution level based on the number of parsed fields.
        - 4: Non-coding variant (e.g., A*01:01:01:02)
        - 3: Synonymous variant (e.g., A*01:01:02)
        - 2: Specific allele (e.g., A*01:01)
        - 1: Allele group (e.g., A*01)
        - 0: Locus only or undefined
        """
        if self.non_coding_variant is not None:
            return 4
        if self.synonymous_variant is not None:
            return 3
        if self.allele is not None:
            return 2
        if self.allele_group is not None:
            return 1
        return 0

    def __eq__(self, other):
        if not isinstance(other, HLA):
            return NotImplemented
        return (
            self.locus == other.locus and
            self.drb_sub_locus == other.drb_sub_locus and
            self.allele_group == other.allele_group and
            self.allele == other.allele and
            self.synonymous_variant == other.synonymous_variant and
            self.non_coding_variant == other.non_coding_variant and
            self.suffix == other.suffix and
            self.group_code == other.group_code
            # TODO: or (self.reduced = other.reduced )
        )

    def __hash__(self):
        return hash(
            (
                self.locus,
                self.drb_sub_locus,
                self.allele_group,
                self.allele,
                self.synonymous_variant,
                self.non_coding_variant,
                self.suffix,
                self.group_code
            )
        )

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
