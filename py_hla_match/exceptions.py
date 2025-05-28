from typing import Union, Dict, Optional


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
            f"You may only compare loci of DRBX. "
            f"Potential error in data preprocessing."
        )
        super().__init__(message)
        self.locus1 = locus1
        self.locus2 = locus2

    def __str__(self):
        return super().__str__()


class MalformedHLADataSourceError(Exception):
    """
    Error raised when an HLA data source is malformed or cannot be parsed.
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


class DataLoaderError(Exception):
    """Base class for exceptions in the loader module."""
    pass


class FileNotFoundError(DataLoaderError, FileNotFoundError):
    """Custom exception for file not found."""
    pass


class UnsupportedFileTypeError(DataLoaderError):
    """Exception for unsupported file types."""
    pass


class EmptyDataError(DataLoaderError):
    """Exception for empty or unparsable data files."""
    pass


class ParsingError:
    """Container for parsing error information."""

    def __init__(
        self,
        row_id: Union[str, int],
        error_type: str,
        message: str,
        details: Optional[Dict] = None
    ):
        self.row_id = row_id
        self.error_type = error_type
        self.message = message
        self.details = details or {}

    def __str__(self):
        return (
            f"Row {self.row_id}: {self.error_type} - {self.message}"
        )
