import logging
import re

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
        # TODO: create private function to check for general syntax errors, e.g. missing seperator
        # -> create a regex that checks for legal HLA strings
        # -> create MalformedHLAStringException
        try:
            # in case of 'HLA-' prefix -> Remove
            allele_string_no_prefix = re.sub(r"HLA-", "", self.allele_string, count=1)
            # split at seperator '*' to optain gene/locus
            allele_string_split_gene = allele_string_no_prefix.split('*')
            self.locus = allele_string_split_gene[0]
            fields_string = allele_string_split_gene[1]
            field_contents = fields_string.split[':']

        except Exception as e:
            # TODO: Add specific Exception
            logger.error(f"Failed to parse HLA allele '{self.allele_str}': {e}")

    def _validate_nomenclature(self) -> None:   
        # TODO: implement
        pass

    def __repr__(self):
        return (
            f"HLA(allele_str='{self.allele_str}', locus='{self.locus}', "
            f"allele_group='{self.allele_group}', allele='{self.allele}', "
            f"synonymous_variant='{self.synonymous_variant}', "
            f"non_coding_variant='{self.non_coding_variant}', suffix='{self.suffix}')"
        )