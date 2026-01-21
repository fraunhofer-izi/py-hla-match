# models.py
import logging
from typing import List
from collections import Counter

from py_hla_match.exceptions import InvalidLocusComparisonError
from py_hla_match.hla import HLA

logger = logging.getLogger(__name__)


class HLAPair:
    def __init__(self, hla1: HLA, hla2: HLA) -> None:
        """
        Pair of HLA objects from the same locus.

        :param hla1: The first HLA object
        :param hla2: The second HLA object
        :raises TypeError: If hla1 or hla2 are not HLA objects
        :raises LocusMismatchError: If HLA objects do not share locus
        """

        # check for object validity -> hla1 and hla2 must be hla objects
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

        self.hla1 = hla1
        self.hla2 = hla2
        self.locus: str = hla1.locus

    @property
    def alleles(self) -> tuple[HLA, HLA]:
        """Returns HLA objects in HLAPair."""
        return (self.hla1, self.hla2)

    # helper functions
    def get_paired_resolution(self) -> int:
        """
        Resolution of HLAPair.
        """
        return min(
            self.hla1.has_resolution_level(),
            self.hla2.has_resolution_level()
        )

    def __str__(self) -> str:
        """String representation of HLAPair."""
        if self.locus is None:
            return "HLAPair(locus=None, hla1=None, hla2=None)"
        hla1_str = str(self.hla1.allele_string)
        hla2_str = str(self.hla2.allele_string)
        return f"HLAPair(locus={self.locus}, hla1={hla1_str}, hla2={hla2_str})"

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other) -> bool:
        if not isinstance(other, HLAPair):
            return NotImplemented
        return frozenset(self.alleles) == frozenset(other.alleles)

    def __hash__(self) -> int:
        return hash(frozenset(self.alleles))


class Individual:
    def __init__(self, hla_data: List[HLAPair]) -> None:
        """
        Represents an individual with HLA data.

        :param hla_data: List of HLAPair
        """
        for item in hla_data:
            if not isinstance(item, HLAPair):
                raise TypeError(
                    f"hla_data list must contain only HLAPair objects, "
                    f"but found an item of type {type(item).__name__}."
                )
        self.hla_data = hla_data
        self._sanity_check()

    def _sanity_check(self) -> None:
        """
        Verify each HLA locus appears only once per individual
        """
        loci_counts = Counter([hla_pair.locus for hla_pair in self.hla_data])
        duplicate_loci = [
            locus for locus, count in loci_counts.items() if count > 1
        ]
        # raise error and report duplicate locus if present
        if duplicate_loci:
            raise ValueError(
                f"Multiple {duplicate_loci} found. Individuals may not "
                f"have multiple HLAPair for the same locus."
            )

    def get_hla_summary(self) -> dict:
        """Get summary of HLA data."""
        if not self.hla_data:
            return {
                "total_loci_typed": 0,
                "resolution_summary": {}
            }

        resolution_levels = [
            pair.get_paired_resolution() for pair in self.hla_data
        ]
        resolution_counts = Counter(resolution_levels)

        return {
            "total_loci_typed": len(self.hla_data),
            "resolution_summary": dict(resolution_counts)
        }


class Patient(Individual):
    """
    Represents an individual in the 'patient' role in research datasets.
    Inherits from Individual.
    """

    def __init__(self,  hla_data: list[HLAPair]) -> None:
        super().__init__(hla_data)


class Donor(Individual):
    """
    Represents an individual in the 'donor' role in research datasets.
    Inherits from Individual.
    """

    def __init__(self,  hla_data: list[HLAPair]) -> None:
        super().__init__(hla_data)
