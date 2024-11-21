import logging
import re

from py_hla_match.exceptions import MalformedHLAStringError

logger = logging.getLogger(__name__)

class HLA:
    def __init__(self, allele_string: str) -> None:
        """
        Initializes an HLA object by parsing the HLA allele string.

        :param allele_string: The HLA allele string to be parsed (e.g., 'A*32:11Q').
        """
        self.allele_string = allele_string
        self.locus = None
        self.allele_group = None
        self.allele = None
        self.synonymous_variant = None
        self.non_coding_variant = None
        self.suffix = None

        self._parse_allele()

    def _parse_allele(self) -> None:
        """
        Parses the HLA allele string and populates the attributes.
        """
        self._validate_nomenclature()

        # Remove 'HLA-' prefix if present
        allele_string_no_prefix = re.sub(r"^HLA-", "", self.allele_string, count=1)

        # Split at '*' to separate locus and fields
        allele_string_split = allele_string_no_prefix.split('*')
        self.locus = allele_string_split[0]
        fields_string = allele_string_split[1]

        # Extract suffix if present (e.g., 'Q', 'N', 'L')
        suffix_match = re.search(r'[A-Z]$', fields_string)
        if suffix_match:
            self.suffix = suffix_match.group(0)
            fields_string = fields_string[:-1]

        # Split fields string at ':' to extract allele components
        field_contents = fields_string.split(':')
        if len(field_contents) > 0:
            self.allele_group = field_contents[0]
        if len(field_contents) > 1:
            self.allele = field_contents[1]
        if len(field_contents) > 2:
            self.synonymous_variant = field_contents[2]
        if len(field_contents) > 3:
            self.non_coding_variant = field_contents[3]

    def _validate_nomenclature(self) -> None:
        """
        Validates the HLA allele string format.
        """
        pattern = r"^(HLA-)?[A-Z0-9]+[*]\d+(:\d+)*[A-Z]?$"
        if not re.match(pattern, self.allele_string):
            raise MalformedHLAStringError(f"Invalid HLA allele string: {self.allele_string}")

    def __repr__(self):
        return (
            f"HLA(allele_string='{self.allele_string}', locus='{self.locus}', "
            f"allele_group='{self.allele_group}', allele='{self.allele}', "
            f"synonymous_variant={self.synonymous_variant}, "
            f"non_coding_variant={self.non_coding_variant}, suffix='{self.suffix}')"
        )
