# parser.py

import logging
from typing import Tuple, List, Dict, Optional, Union, Any
import pandas as pd

from py_hla_match.hla import HLA
from py_hla_match.models import Individual, HLAPair
from py_hla_match.exceptions import (
    MalformedHLAStringError,
    MalformedHLADataSourceError,
    ParsingError
)


logger = logging.getLogger(__name__)


class HLAParser:
    """
    Parses HLA data from HLADataSource.

    Parses  HLA allele strings from DataFrame(s) into interna Individual
    objects.

    CONVENTIONS:
    - single row = single individual's HLA data

    Supports three data structures:
    - 'single': One individual per row
    - 'paired': Patient and donor in same row
    - 'panel': Multiple individuals (e.g., family) in same row
    """

    def __init__(self):
        """Initialize parser with error collection."""
        self.errors: List[ParsingError] = []
        self._reset_errors()

    def parse(
        self,
        df: pd.DataFrame,
        structure: str,
        column_mapping: Dict[str, List[str]],
        id_column: Optional[str] = None
    ) -> Union[
        # ID and single individuals
        List[Tuple[Union[str, int], Individual]],
        # ID and patient/donor pairs
        List[Tuple[Union[str, int], Individual, Individual]],
        # ID and panel members
        List[Tuple[Union[str, int], Dict[str, Individual]]]
    ]:
        """
        Parse DataFrame into Individual objects based on structure.

        Args:
            df: DataFrame containing HLA data
            structure: Data structure - 'single', 'paired', or 'panel'
            column_mapping: Dict mapping to HLA columns
                - For 'single': {'individual': [columns]}
                - For 'paired': {'patient': [cols], 'donor': [cols]}
                - For 'panel': {'member1': [cols], 'member2': [cols], ...}
            id_column: Column with IDs (None = use index/'_row_index')

        Returns:
            Based on structure:
            - 'single': List[(id, Individual)]
            - 'paired': List[(id, Patient, Donor)]
            - 'panel': List[(id, {'name': Individual, ...})]

        Raises:
            ValueError: Invalid structure or column_mapping
            MalformedHLADataSourceError: Critical data structure issues

        Example:
            >>> parser = HLAParser()
            >>> pairs = parser.parse(
            ...     df,
            ...     structure='paired',
            ...     column_mapping={
            ...         'patient': ['pat_A1', 'pat_A2'],
            ...         'donor': ['don_A1', 'don_A2']
            ...     }
            ... )
        """
        # reset erros
        self._reset_errors()

        # validate
        self._validate_inputs(structure, column_mapping, df)

        # get ID if ID
        if id_column and id_column not in df.columns:
            raise ValueError(
                f"ID column '{id_column}' not found in DataFrame. "
                f"Available columns: {list(df.columns)}"
            )

        # select parser for structure
        if structure == 'single':
            return self._parse_single(df, column_mapping, id_column)
        elif structure == 'paired':
            return self._parse_paired(df, column_mapping, id_column)
        elif structure == 'panel':
            return self._parse_panel(df, column_mapping, id_column)
        else:
            raise ValueError(
                f"Invalid structure '{structure}'. "
                "Must be 'single', 'paired', or 'panel'"
            )

    def _create_hla_pairs_from_strings(
            self, hla_strings: List[str], context: Dict[str, Any]
    ) -> List[HLAPair]:
        """
        HLA strings to HLAPairs grouped by locus.

        Args:
            hla_strings: List of HLA allele strings from row
            context: Error context

        Returns:
            List of HLAPair objects
        """
        # TODO: Placeholder for future implementation
        pass

    def get_error_summary(self) -> str:
        """Summary of parsing errors."""
        if not self.errors:
            return "No parsing errors encountered."

        error_types = {}
        for error in self.errors:
            error_types[error.error_type] = (
                error_types.get(error.error_type, 0) + 1
            )

        summary = f"Total errors: {len(self.errors)}\n"
        for error_type, count in error_types.items():
            summary += f"  - {error_type}: {count}\n"

        return summary

    def _reset_errors(self) -> None:
        """Resets error collection for new session."""
        self.errors = []

    def _validate_inputs(
        self,
        structure: str,
        column_mapping: Dict,
        df: pd.DataFrame
    ) -> None:
        """Validate parsing inputs."""
        # valid structures
        valid_structures = {'single', 'paired', 'panel'}
        if structure not in valid_structures:
            raise ValueError(
                f"Invalid structure '{structure}'. "
                f"Must be one of: {valid_structures}"
            )

        # validate column mapping structure
        if structure == 'single' and 'individual' not in column_mapping:
            raise ValueError(
                "For 'single' structure, column_mapping must "
                "contain 'individual' key"
            )
        elif structure == 'paired':
            if not all(
                key in column_mapping for key in ['patient', 'donor']
            ):
                raise ValueError(
                    "For 'paired' structure, column_mapping must "
                    "contain 'patient' and 'donor' keys"
                )

        # Check that specified columns exist
        all_columns = []
        for columns in column_mapping.values():
            all_columns.extend(columns)

        missing_columns = set(all_columns) - set(df.columns)
        if missing_columns:
            raise ValueError(
                f"Columns not found in DataFrame: {missing_columns}"
            )

    def _get_row_id(
        self,
        row: pd.Series,
        idx: int,
        id_column: Optional[str]
    ) -> Union[str, int]:
        """ID from specified column else index."""
        if id_column:
            return row[id_column]
        elif '_row_index' in row:
            return row['_row_index']
        else:
            return idx

    def _parse_single(
        self,
        df: pd.DataFrame,
        column_mapping: Dict[str, List[str]],
        id_column: Optional[str]
    ) -> List[Tuple[Union[str, int], Individual]]:
        """Parse single individual."""
        results = []
        hla_columns = column_mapping['individual']

        for idx, row in df.iterrows():
            row_id = self._get_row_id(row, idx, id_column)

            try:
                individual = self._parse_individual_from_row(
                    row[hla_columns],
                    row_id,
                    'individual'
                )
                results.append((row_id, individual))

            except Exception as e:
                self._record_error(
                    row_id,
                    'ParseError',
                    f"Failed to parse individual: {e}"
                )
                # continue processing rows

        self._log_parsing_summary(len(df), len(results))
        return results

    def _parse_paired(
        self,
        df: pd.DataFrame,
        column_mapping: Dict[str, List[str]],
        id_column: Optional[str]
    ) -> List[Tuple[Union[str, int], Individual, Individual]]:
        """Parse patient/donor pairs."""
        results = []
        patient_columns = column_mapping['patient']
        donor_columns = column_mapping['donor']

        for idx, row in df.iterrows():
            row_id = self._get_row_id(row, idx, id_column)

            try:
                patient = self._parse_individual_from_row(
                    row[patient_columns],
                    row_id,
                    'patient'
                )
                donor = self._parse_individual_from_row(
                    row[donor_columns],
                    row_id,
                    'donor'
                )
                results.append((row_id, patient, donor))

            except Exception as e:
                self._record_error(
                    row_id,
                    'ParseError',
                    f"Failed to parse pair: {e}"
                )
                # continue processing rows

        self._log_parsing_summary(len(df), len(results))
        return results

    def _parse_panel(
        self,
        df: pd.DataFrame,
        column_mapping: Dict[str, List[str]],
        id_column: Optional[str]
    ) -> List[Tuple[Union[str, int], Dict[str, Individual]]]:
        """Parse multiple individuals (panel)."""
        results = []

        for idx, row in df.iterrows():
            row_id = self._get_row_id(row, idx, id_column)
            panel = {}

            try:
                for member_name, columns in column_mapping.items():
                    individual = self._parse_individual_from_row(
                        row[columns],
                        row_id,
                        member_name
                    )
                    panel[member_name] = individual

                results.append((row_id, panel))

            except Exception as e:
                self._record_error(
                    row_id,
                    'ParseError',
                    f"Failed to parse panel: {e}"
                )
                # continue processing rows

        self._log_parsing_summary(len(df), len(results))
        return results

    # the actual parsing logic
    def _parse_individual_from_row(
        self,
        hla_values: pd.Series,
        row_id: Union[str, int],
        individual_type: str
    ) -> Individual:
        """Parse HLA values into an Individual object."""
        hla_pairs = []
        locus_alleles = {}

        # Parse each HLA string
        for column_name, hla_string in hla_values.items():
            # Skip empty values
            if not hla_string or str(hla_string).strip() == '':
                continue

            try:
                # HLA object
                hla = HLA(hla_string)

                # group by locus
                if hla.locus not in locus_alleles:
                    locus_alleles[hla.locus] = []
                locus_alleles[hla.locus].append(hla)

            except (MalformedHLAStringError, Exception):
                self._record_error(
                    row_id,
                    'InvalidHLAString',
                    f"Invalid HLA string '{hla_string}' in column "
                    f"'{column_name}' for {individual_type}",
                    details={
                        'column': column_name,
                        'value': hla_string,
                        'individual_type': individual_type
                    }
                )

        # create HLAPair objects for locus
        for locus, alleles in locus_alleles.items():
            if len(alleles) == 2:
                hla_pairs.append(HLAPair(alleles[0], alleles[1]))
            elif len(alleles) == 1:
                # Single allele at locus - create pair with None
                hla_pairs.append(HLAPair(alleles[0], None))
                self._record_error(
                    row_id,
                    'IncompleteLocus',
                    f"Only one allele found for locus {locus} "
                    f"in {individual_type}",
                    details={
                        'locus': locus,
                        'individual_type': individual_type
                    }
                )
            elif len(alleles) > 2:
                # Too many alleles - data structure issue
                raise MalformedHLADataSourceError(
                    f"Found {len(alleles)} alleles for locus {locus} "
                    f"in {individual_type} at row {row_id}. "
                    f"Maximum 2 alleles per locus allowed."
                )
            # TODO: do we need to handle cases with no alleles?

        return Individual(hla_pairs)

    def _record_error(
        self,
        row_id: Union[str, int],
        error_type: str,
        message: str,
        details: Optional[Dict] = None
    ) -> None:
        """Record a parsing error."""
        error = ParsingError(row_id, error_type, message, details)
        self.errors.append(error)
        logger.warning(str(error))

    def _log_parsing_summary(
        self,
        total_rows: int,
        successful_rows: int
    ) -> None:
        """Log summary of parsing results."""
        failed_rows = total_rows - successful_rows

        logger.info(
            f"Parsing complete: {successful_rows}/{total_rows} rows "
            f"successfully parsed"
        )
        if failed_rows > 0:
            logger.warning(
                f"{failed_rows} rows failed to parse. "
                "Check logs for details."
            )
        if self.errors:
            logger.warning(
                f"Encountered {len(self.errors)} parsing errors. "
                "Call get_error_summary() for details."
            )
