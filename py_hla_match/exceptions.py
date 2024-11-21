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