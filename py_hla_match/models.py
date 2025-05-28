# models.py
import logging
from typing import List, Optional
from collections import Counter

from py_hla_match.exceptions import InvalidLocusComparisonError
from py_hla_match.hla import HLA

logger = logging.getLogger(__name__)


class HLAPair:
    def __init__(self, hla1: Optional[HLA], hla2: Optional[HLA]) -> None:
        """
        Pair of HLA objects from the same locus.

        :param hla1: The first HLA object (None if unparsable)
        :param hla2: The second HLA object (None if unparsable)
        """

        # check for object validity -> hla1 and hla2 must be hla objects
        if hla1 is not None and not isinstance(hla1, HLA):
            raise TypeError(
                f"hla1 must be an instance of HLA, not {type(hla1).__name__}."
            )
        if hla2 is not None and not isinstance(hla2, HLA):
            raise TypeError(
                f"hla2 must be an instance of HLA, not {type(hla2).__name__}."
            )

        self.hla1 = hla1
        self.hla2 = hla2
        self.locus = self._get_locus()

    def _get_locus(self) -> Optional[str]:
        """
        Returns the locus of the HLA pair.
        """
        locus1 = self.hla1.locus if self.hla1 else None
        locus2 = self.hla2.locus if self.hla2 else None

        # if both are None, return None
        if locus1 is None and locus2 is None:
            return None

        # handle DRBX fist
        # may be obsolete since we adjusted hla.py
        locus_drbx = [
            'DRB2', 'DRB3', 'DRB4', 'DRB5', 'DRB6', 'DRB7', 'DRB8', 'DRB9',
            'DRBX'
        ]

        # for now just hardcoded DRBX
        # TODO: needs external validation
        if locus1 in locus_drbx:
            return 'DRBX'
        elif locus2 in locus_drbx:
            return 'DRBX'

        # if not drbx HLAPair should be considered from the same locus
        if locus1 is None:
            return locus2
        if locus2 is None:
            return locus1

        if locus1 == locus2:
            return locus1
        else:
            raise InvalidLocusComparisonError(locus1, locus2)

    # helper functions
    def has_any_data(self) -> bool:
        """Check if HLAPair has any HLA data."""
        return self.hla1 is not None or self.hla2 is not None

    def has_locus_info(self) -> bool:
        """Check if locus is available."""
        return self.locus is not None

    def has_valid_hla(self) -> bool:
        """
        Check if at least one allele has valid HLA data (one-field or better).
        """
        hla1_valid = self.hla1 and self.hla1.has_one_field_resolution()
        hla2_valid = self.hla2 and self.hla2.has_one_field_resolution()
        return hla1_valid or hla2_valid

    def has_high_resolution_hla(self) -> bool:
        """
        Check if at least one allele has high-resolution HLA data (two-field or
        better).
        """
        hla1_high_res = self.hla1 and self.hla1.has_two_field_resolution()
        hla2_high_res = self.hla2 and self.hla2.has_two_field_resolution()
        return hla1_high_res or hla2_high_res

    def has_hla_pair(self) -> bool:
        """Check if both alleles are present."""
        return self.hla1 is not None and self.hla2 is not None

    def has_valid_hla_pair(self) -> bool:
        """Check if both alleles have valid HLA data (one-field or better)."""
        hla1_valid = self.hla1 and self.hla1.has_one_field_resolution()
        hla2_valid = self.hla2 and self.hla2.has_one_field_resolution()
        return hla1_valid and hla2_valid

    def has_high_resolution_hla_pair(self) -> bool:
        """
        Check if both alleles have high-resolution HLA data (two-field or
        better).
        """
        hla1_high_res = self.hla1 and self.hla1.has_two_field_resolution()
        hla2_high_res = self.hla2 and self.hla2.has_two_field_resolution()
        return hla1_high_res and hla2_high_res

    def __str__(self) -> str:
        """String representation of HLAPair."""
        if self.locus is None:
            return "HLAPair(locus=None, hla1=None, hla2=None)"
        hla1_str = str(self.hla1.allele_string) if self.hla1 else "None"
        hla2_str = str(self.hla2.allele_string) if self.hla2 else "None"
        return f"HLAPair(locus={self.locus}, hla1={hla1_str}, hla2={hla2_str})"

    def __repr__(self) -> str:
        return self.__str__()


class Individual:
    def __init__(self, hla_data: List[HLAPair]) -> None:
        """
        Represents an individual with HLA data.

        :param hla_data: List of HLAPair objects
        """
        self.hla_data = hla_data
        self._sanity_check()

    def _sanity_check(self) -> None:
        """
        Verify each HLA locus appreas only once in each individual
        """
        valid_loci = [pair for pair in self.hla_data if pair.locus is not None]

        loci_counts = Counter([hla_pair.locus for hla_pair in valid_loci])
        duplicate_loci = [
            locus for locus, count in loci_counts.items() if count > 1
        ]
        # raise error and report duplicate locus if present
        if duplicate_loci:
            raise ValueError(
                f"Duplicate loci found: {duplicate_loci}. Individuals may not "
                f"have mutliple HLA pairs for the same locus."
            )

    def get_hla_summary(self) -> dict:
        """Get summary of HLA data."""
        total_pairs = len(self.hla_data)
        parsed_pairs = sum(
            1 for pair in self.hla_data if pair.has_any_data()
        )
        valid_pairs = sum(
            1 for pair in self.hla_data if pair.has_valid_hla_pair()
        )
        high_res_pairs = sum(
            1 for pair in self.hla_data if pair.has_high_resolution_hla_pair()
        )

        return {
            "total_alles": total_pairs,
            "parsed_alleles": parsed_pairs,
            "valid_alleles": valid_pairs,
            "high_resolution_alleles": high_res_pairs
        }

    def get_best_match(self, individuals: List['Individual']) -> 'Donor':
        """
        Get the best match from a list of donors.
        """
        pass


class Patient(Individual):
    """
    Represents a patient, inheriting from Individual.
    """

    def __init__(self,  hla_data: list[HLAPair]) -> None:
        super().__init__(hla_data)


class Donor(Individual):
    """
    Represents a donor, inheriting from Individual.
    """

    def __init__(self,  hla_data: list[HLAPair]) -> None:
        super().__init__(hla_data)
