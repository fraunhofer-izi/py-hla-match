import logging
from .hla import HLA

logger = logging.getLogger(__name__)

class MatchResult():
    """
    Report class on patient/donor compatability.
    """

    def __init__(self, match_type: str, homozygous_patient: bool) -> None:
        self.match_type = match_type
        self.homozygous_patient = homozygous_patient
        self.details = self._get_details()
        self.match_score = self._compute_score()

    def _get_details(self) -> str:
        if self.match_type == "full_match":
            return "Three-field match: Match including the synonymous mutation."
        if self.match_type == "pot_mutation_mismatch":
            return "At least two-field match: Match including the specific HLA allele, no information about the synonymous mutation."
        if self.match_type == "mutation_mismatch":
            return "Two-field match: Match including the specific HLA allele, synonymous mutation."
        if self.match_type == "pot_allele_mismatch":
            return "At least one-field match: Match, no information about the specific HLA allele."
        if self.match_type == "allele_mismatch":
            return "Two-field mismatch: Mismatch of the specific HLA allele."
        if self.match_type == "antigen_mismatch":
            return "One-field mismatch: Mismatch of the encoded antigen."
        return "Unknown match type"


    def _compute_score(self):
        """
        TODO: base on clinicians feedback
        """
        return None
    

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
        
        if patient.hla1.allele_string == patient.hla2.allele_string:
            homozygous_patient = True
        else:
            homozygous_patient = False

        # TODO: implement
        allele_group_match = None
        allele_match = None
        mutation_group_match = None

        if allele_group_match and allele_match and mutation_group_match:
            match_type = "full_match"
        elif allele_group_match and allele_match and (mutation_group_match is None):
            match_type = "potential_mutation_mismatch"
        elif allele_group_match and allele_match and not mutation_group_match:
            match_type = "mutation_mismatch"
        elif allele_group_match and (allele_match is None):
            match_type = "pot_allele_mismatch"
        elif allele_group_match and not allele_match:
            match_type = "allele_mismatch"
        elif not allele_group_match:
            match_type = "antigen_mismatch"
        else:
            match_type = "NA"

        match_result = MatchResult(match_type=match_type, homozygous_patient=homozygous_patient)

        return match_result
        


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


