import logging
import re
import pyard

from py_hla_match.exceptions import MalformedHLAStringError

logger = logging.getLogger(__name__)


class HLA:
    def __init__(self, allele_string: str) -> None:
        """
        Initializes an HLA object by parsing the HLA allele string.

        :param allele_string: The HLA allele string to be parsed
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

        self._parse_allele()
        self._ard_redux()

    def _parse_allele(self) -> None:
        """
        Parses the HLA allele string and populates the attributes.
        """
        # validate the allele string
        match = self._validate_nomenclature()

        # extract locus, allele fields, suffix or group code
        self.locus = match.group('locus')
        allele_fields = match.group('allele_fields')
        self.suffix = match.group('suffix')
        self.group_code = match.group('group_code')

        # extract detailes from allele fields
        field_contents = allele_fields.split(':')
        if len(field_contents) > 0:
            self.allele_group = field_contents[0]
        if len(field_contents) > 1:
            self.allele = field_contents[1]
        if len(field_contents) > 2:
            self.synonymous_variant = field_contents[2]
        if len(field_contents) > 3:
            self.non_coding_variant = field_contents[3]

    def _validate_nomenclature(self) -> re.Match:
        """
        Validates the HLA allele string format.
        """
        pattern = re.compile(
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
        match = pattern.match(self.allele_string)
        if not match:
            raise MalformedHLAStringError(
                f"Invalid HLA allele string: {self.allele_string}"
            )
        # TODO: Return additional details on specific error occurence
        # TODO: verify locus is part of named genes within the HLA region
        # TODO: potentially validate the allele string is known allele with
        # ref. to https://hla.alleles.org/alleles/index.html)
        return match

    def _ard_redux(self):
        """
        Reduce Allele to 2 field ARD level using the wonderful py-ard package
        For further information: https://github.com/nmdp-bioinformatics/py-ard
        """
        ard = pyard.init()
        redux_type = 'lgx'
        self.ard_redux_allele_string = ard.redux(
            self.allele_string, redux_type
        )
        pattern = re.compile(
            r"""
            ^(?:HLA-)?                 # Optional 'HLA-' prefix
            (?P<locus>[A-Z0-9]+)       # Locus
            \*                         # Asterisk separator
            (?P<allele_fields>\d{2,}(?::\d{2,}){0,3})  # Allele fields
            """,
            re.VERBOSE
        )
        match = pattern.match(self.ard_redux_allele_string)
        if match:
            allele_fields = match.group('allele_fields')
            field_contents = allele_fields.split(':')
            if len(field_contents) > 0:
                self.ard_redux_allele_group = field_contents[0]
            if len(field_contents) > 1:
                self.ard_redux_allele = field_contents[1]

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


if __name__ == "__main__":
    hla = HLA('A*02:01:02')
    print(hla)
