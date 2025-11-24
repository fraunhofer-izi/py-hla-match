import logging
import re
import threading
from typing import Optional

from pyard.exceptions import InvalidAlleleError

from py_hla_match.singleton import get_ard_instance
from py_hla_match.config import (
    get_config,
    get_config_version
)
from py_hla_match.exceptions import (
    MalformedHLAStringError,
    EmptyHLAStringError, PyardLibraryError
)

logger = logging.getLogger(__name__)

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


class HLA:
    _cache = {}
    _redux_cache = {}
    _cache_lock = threading.RLock()
    _redux_cache_lock = threading.RLock()
    _config_version: Optional[int] = None

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

    def __new__(cls, allele_string: str):
        """Thread-safe caching"""
        with cls._cache_lock:
            # Invalidate cache if config version changed
            current_ver = get_config_version()
            if cls._config_version != current_ver:
                cls._cache.clear()
                cls._redux_cache.clear()
                cls._config_version = current_ver
            if allele_string in cls._cache:
                return cls._cache[allele_string]

            # Create new instance only if not cached
            instance = super().__new__(cls)
            cls._cache[allele_string] = instance
            return instance

    def __init__(self, allele_string: str) -> None:
        # Skip initialization if already initialized (cached object)
        if hasattr(self, '_locked') and self._locked:
            return
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

        # extract locus
        self.locus = match.group('locus')
        # handle DRB3/4/5 region
        if self.locus in get_config().drb345_sub_loci:
            # if locus is DRB3/4/5, set the sub-locus
            self.drb_sub_locus = self.locus
            self.locus = 'DRB345'

        allele_fields = match.group('allele_fields')
        nan_field = match.group('nan')
        remainder = match.group('remainder')

        if allele_fields:
            self.suffix = match.group('suffix')
            self.group_code = match.group('group_code')
            if self.group_code == 'G' and allele_fields.count(':') != 2:
                raise MalformedHLAStringError(
                    f"'{self.allele_string}' – 'G' group must have 3 fields."
                )
            if self.group_code == 'P' and allele_fields.count(':') != 1:
                raise MalformedHLAStringError(
                    f"'{self.allele_string}' – 'P' group must have 2 fields."
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
            raise EmptyHLAStringError(
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

        match = get_config().nomenclature_pattern.match(self.allele_string)
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
        return locus in get_config().effective_valid_loci

    def _ard_redux(self):
        """Thread-safe ARD redux with caching"""
        redux_type = 'lgx'
        cache_key = (self.allele_string, redux_type)

        # Thread-safe cache check
        with self._redux_cache_lock:
            if cache_key in self._redux_cache:
                cached_result = self._redux_cache[cache_key]
                self.ard_redux_allele_string = cached_result['redux_string']
                self.ard_redux_allele_group = cached_result['allele_group']
                self.ard_redux_allele = cached_result['allele']
                return

        try:
            ard = get_ard_instance()
            redux_string = ard.redux(self.allele_string, redux_type).strip()
        except InvalidAlleleError as e:
            # propagate allele specific exception
            raise e
        except Exception as e:
            # catch and re-raise any other (library specific) exceptions
            raise PyardLibraryError(
                f"Failed during allele reduction for '{self.allele_string}' with redux_type '{redux_type}'.",
                details=str(e)
            ) from e

        # Parse the result
        match = REDUX_PATTERN.match(redux_string)
        if not match:
            raise MalformedHLAStringError(
                "py-ard reports unexpected string not matching regex "
                f"'{redux_string}'."
            )

        allele_group = None
        allele = None
        if match:
            allele_fields = match.group('allele_fields')
            field_contents = allele_fields.split(':')
            if len(field_contents) > 0:
                allele_group = field_contents[0]
            if len(field_contents) > 1:
                allele = field_contents[1]

        # Thread-safe caching
        with self._redux_cache_lock:
            # Double-check pattern: another thread might have computed this
            # while we were working
            if cache_key in self._redux_cache:
                cached_result = self._redux_cache[cache_key]
                self.ard_redux_allele_string = cached_result['redux_string']
                self.ard_redux_allele_group = cached_result['allele_group']
                self.ard_redux_allele = cached_result['allele']
            else:
                # We're the first to compute this, cache it
                self.ard_redux_allele_string = redux_string
                self.ard_redux_allele_group = allele_group
                self.ard_redux_allele = allele

                self._redux_cache[cache_key] = {
                    'redux_string': redux_string,
                    'allele_group': allele_group,
                    'allele': allele
                }

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
