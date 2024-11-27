from py_hla_match.models import Patient, Donor


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
    
def match(patient: Patient, donor: Donor) -> MatchResult:

        if patient.hla1 == patient.hla2:
            homozygous_patient = True
        else:
            homozygous_patient = False

        # TODO: implement
        allele_group_match = None
        allele_match = None
        mutation_group_match = None

        if not allele_group_match:
            match_type = "antigen_mismatch"
        elif allele_group_match and not allele_match:
            match_type = "allele_mismatch"
        elif allele_group_match and (allele_match is None):
            match_type = "pot_allele_mismatch"
        elif allele_group_match and allele_match and mutation_group_match:
            match_type = "full_match"
        elif allele_group_match and allele_match and not mutation_group_match:
            match_type = "mutation_mismatch"
        elif allele_group_match and allele_match and (mutation_group_match is None):
            match_type = "potential_mutation_mismatch"
        else:
            match_type = "NA"

        match_result = MatchResult(match_type=match_type, homozygous_patient=homozygous_patient)

        return match_result