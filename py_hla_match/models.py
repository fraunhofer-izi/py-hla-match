import logging

from py_hla_match.hla import HLA

logger = logging.getLogger(__name__)


class HLAPair:
    def __init__(self, hla1: HLA, hla2: HLA, locus: str) -> None:
        """
        Initializes an HLAPair object with two HLA objects and a locus.

        :param hla1: The first HLA object
        :param hla2: The second HLA object
        :param locus: The locus of the HLA pair
        """

        self.hla1 = hla1
        self.hla2 = hla2
        self.locus = locus


class Individual:
    def __init__(self, hla1: HLA, hla2: HLA) -> None:
        """
        Initializes an individual with two HLA objects.

        :param hla1: The first HLA object
        :param hla2: The second HLA object
        """

        self.hla1 = hla1
        self.hla2 = hla2


class Patient(Individual):
    """
    Represents a patient, inheriting from Individual.
    """

    def __init__(self, hla1: HLA, hla2: HLA) -> None:
        super().__init__(hla1, hla2)


class Donor(Individual):
    """
    Represents a donor, inheriting from Individual.
    """

    def __init__(self, hla1: HLA, hla2: HLA) -> None:
        super().__init__(hla1, hla2)
