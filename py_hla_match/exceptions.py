class MalformedHLAStringError(Exception):
    """
    Error raised when an HLA string is malformed or cannot be parsed.
    """
    def __init__(self, message, details=None):
        """
        :param message: Error message.
        :param details: Optional additional details about the error.
        """
        super().__init__(message)
        self.details = details

    def __str__(self):
        base_message = super().__str__()
        if self.details:
            return f"{base_message} (Details: {self.details})"
        return base_message


class InvalidLocusComparisonError(Exception):
    """
    Exception raised when attempting to compare alleles from different loci.
    """
    def __init__(self, locus1, locus2):
        message = (
            f"Invalid locus comparison between '{locus1}' and '{locus2}'. "
            f"You may only compare DRBX. "
            f"Potential error in data preprocessing."
        )
        super().__init__(message)
        self.locus1 = locus1
        self.locus2 = locus2

    def __str__(self):
        return super().__str__()
