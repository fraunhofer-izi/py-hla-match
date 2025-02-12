import logging

from py_hla_match.exceptions import InvalidLocusComparisonError
from py_hla_match.hla import HLA

logger = logging.getLogger(__name__)


class HLAPair:
    def __init__(self, hla1: HLA, hla2: HLA) -> None:
        """
        Pair of HLA objects for the same locus.

        :param hla1: The first HLA object
        :param hla2: The second HLA object
        """

        self.hla1 = hla1
        self.hla2 = hla2
        self.locus = self._get_locus()

    def _get_locus(self) -> str:
        """
        Returns the locus of the HLA pair.
        """
        if self.hla1 is None or self.hla2 is None:
            return None
        elif self.hla1.locus is None or self.hla2.locus is None:
            return None
        elif self.hla1.locus != self.hla2.locus:
            raise InvalidLocusComparisonError("HLA alleles must be from the same locus.")
        else:
            return self.hla1.locus


class Individual:
    def __init__(self, hla_data: [HLAPair]) -> None:
        """
        Represents an individual with HLA data.

        :param hla_data: List of HLAPair objects
        """
        self.hla_data = hla_data


class Patient(Individual):
    """
    Represents a patient, inheriting from Individual.
    """

    def __init__(self,  hla_data: [HLAPair]) -> None:
        super().__init__(hla_data)

    def match(self, donor: 'Donor') -> None:
        """
        Match the patient with a donor and get compatibility.
        """
        pass

    def get_best_match(self, donors: ['Donor']) -> 'Donor':
        """
        Get the best match from a list of donors.
        """
        pass


class Donor(Individual):
    """
    Represents a donor, inheriting from Individual.
    """

    def __init__(self,  hla_data: [HLAPair]) -> None:
        super().__init__(hla_data)
