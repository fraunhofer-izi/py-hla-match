"""
HLA Data Loader Module

Provides simple, robust loading of HLA data from CSV and Excel files.
"""

import logging
import pandas as pd
from pathlib import Path
from typing import Optional, List, Union

from py_hla_match.exceptions import (
    FileNotFoundError,
    DataLoaderError,
    EmptyDataError,
    UnsupportedFileTypeError
)

logger = logging.getLogger(__name__)


class HLADataLoader:
    """
    Load HLA data from CSV or Excel files into pandas DataFrames.

    Simple interface that handles file validation, loading, and
    basic data quality checks.
    """

    SUPPORTED_EXTENSIONS = {'.csv', '.xls', '.xlsx'}

    @staticmethod
    def load(
        filepath: Union[str, Path],
        columns: Optional[List[str]] = None,
        id_column: Optional[str] = None,
        sheet_name: Union[str, int] = 0,
        delimiter: str = ','
    ) -> pd.DataFrame:
        """
        Load HLA data from file into a DataFrame.

        Args:
            filepath: Path to CSV or Excel file
            columns: Specific columns to load (None = all columns)
            id_column: Column with ID (None = use row index)
            sheet_name: Excel files, sheet to load (default: first sheet)
            delimiter: CSV files, column delimiter (default: comma)

        Returns:
            pd.DataFrame: HLA data as DataFrame

        Raises:
            FileNotFoundError: File doesn't exist
            UnsupportedFileTypeError: Unsupported file type
            EmptyDataError: No data found in file
            DataLoaderError: Other loading errors

        Example:
            >>> df = HLADataLoader.load(
            ...     "donors.csv",
            ...     columns=['donor_id', 'HLA_A1', 'HLA_A2'],
            ...     id_column='donor_id'
            ... )
        """
        # basic validation of exist(path, file) and non-empty
        filepath = Path(filepath).resolve()
        HLADataLoader._validate_file(filepath)

        # file type validation
        file_extension = filepath.suffix.lower()
        if file_extension not in HLADataLoader.SUPPORTED_EXTENSIONS:
            raise UnsupportedFileTypeError(
                f"Unsupported file type: '{file_extension}'. "
                f"Supported types: {HLADataLoader.SUPPORTED_EXTENSIONS}"
            )

        try:
            # select loader for file type
            # CSV
            if file_extension == '.csv':
                df = HLADataLoader._load_csv(
                    filepath, columns, delimiter
                )
            else:  # Excel
                df = HLADataLoader._load_excel(
                    filepath, columns, sheet_name
                )

            # non-empty check (e.g. header only)
            if df.empty:
                raise EmptyDataError(f"No data found in file: {filepath}")

            # row index if no ID
            if id_column is None:
                df['_row_index'] = range(len(df))
                logger.info(
                    f"No ID column specified. Added '_row_index' column "
                    f"with values 0 to {len(df)-1}"
                )
            elif id_column not in df.columns:
                raise DataLoaderError(
                    f"Specified ID column '{id_column}' not found in file. "
                    f"Available columns: {list(df.columns)}"
                )

            logger.info(
                f"Successfully loaded {len(df)} rows from {filepath.name}"
            )

            return df

        except (EmptyDataError, DataLoaderError, UnsupportedFileTypeError):
            # custom exceptions are already raised
            raise
        except Exception as e:
            # (unknown) other exceptions
            raise DataLoaderError(
                f"Failed to load file '{filepath}': {type(e).__name__} - {e}"
            ) from e

    @staticmethod
    def _validate_file(filepath: Path) -> None:
        """Validate file exists and is accessible."""
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        if not filepath.is_file():
            raise FileNotFoundError(f"Not a file: {filepath}")
        if not filepath.stat().st_size > 0:
            raise EmptyDataError(f"File is empty: {filepath}")

    @staticmethod
    def _load_csv(
        filepath: Path,
        columns: Optional[List[str]],
        delimiter: str
    ) -> pd.DataFrame:
        """Load CSV file with specified parameters."""
        logger.debug(f"Loading CSV file: {filepath}")

        # minimal assumptions
        df = pd.read_csv(
            filepath,
            usecols=columns,
            delimiter=delimiter,
            dtype=str,  # string for HLA data
            na_filter=False,  # no NaN because they can make trouble
            skip_blank_lines=True
        )

        return df

    @staticmethod
    def _load_excel(
        filepath: Path,
        columns: Optional[List[str]],
        sheet_name: Union[str, int]
    ) -> pd.DataFrame:
        """Load Excel file with specified parameters."""
        logger.debug(
            f"Loading Excel file: {filepath}, sheet: {sheet_name}"
        )

        # minimal assumptions
        df = pd.read_excel(
            filepath,
            sheet_name=sheet_name,
            usecols=columns,
            dtype=str,  # string for HLA data
            na_filter=False  # no NaN because they can make troubles
        )

        return df
