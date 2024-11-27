import logging

from py_hla_match.matching import MatchResult, match
from .hla import HLA

logger = logging.getLogger(__name__)

class Individual:
    def __init__(self, hla1: HLA, hla2: HLA) -> None:
        """
        Initializes an individual with two HLA objects.

        :param hla1: The first HLA object
        :param hla2: The second HLA object
        """
        
        self.hla1 = hla1
        self.hla2 = hla2

    def match(self, other: "Individual") -> MatchResult:
        """
        Matches this individual against another and determines compatibility.

        Compatibility is determined by comparing the HLA alleles (hla1 and hla2) 
        of the two individuals. The match type, homozygous status of the patient, 
        and detailed results are encapsulated in a `MatchResult` object.

        :param other: Another Individual (either a Patient or Donor) to compare against.
        :return: A `MatchResult` object that includes:
            - match_type: A string describing the compatibility level (e.g., full_match, allele_mismatch).
            - homozygous_patient: A boolean indicating if the patient is homozygous.
            - details: A human-readable description of the compatibility assessment.
        """        

        # Determine patient/donor roles based on class identity
        if isinstance(self, Patient) and isinstance(other, Donor):
            patient = self
            donor = other
        elif isinstance(self, Donor) and isinstance(other, Patient):
            patient = other
            donor = self
        else:
            raise ValueError("Matching requires one Patient and one Donor.")
        
        return match(patient, donor)
        

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


