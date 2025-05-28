# export.py

import logging
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Union
import pandas as pd

from py_hla_match.exceptions import DataLoaderError, ParsingError

logger = logging.getLogger(__name__)


def export_results(
    results_df: pd.DataFrame,
    output_file: Union[str, Path],
    errors: Optional[List[ParsingError]] = None,
    write_error_log: bool = True
) -> None:
    """
    Export HLA matching results with optional error logging.

    File format is auto-detected from extension (.csv or .xlsx).
    If write_error_log is True, creates ab error log file.
    Args:
        results_df: DataFrame of HLA matching results
        output_file: Output file path (format detected from extension)
        errors: Parsing errors from HLAParser (optional)
        write_error_log: Bool for writing error log file (or not)

    Raises:
        ValueError: Unsupported file format
        DataLoaderError: Export failed

    Example:
        >>> # Basic export
        >>> export_results(results_df, "matches.xlsx")

        >>> # With error logging from parser
        >>> export_results(results_df, "matches.csv", parser.errors)

        >>> # Without error log
        >>> export_results(results_df, "output.xlsx", write_error_log=False)
    """
    output_path = Path(output_file).resolve()

    # create directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # check if DataFrame is empty
    if results_df.empty:
        logger.warning(f"Exporting empty DataFrame to {output_path}")

    # Export based on file extension
    file_extension = output_path.suffix.lower()

    try:
        if file_extension == '.csv':
            results_df.to_csv(output_path, index=False)
            logger.info(
                f"Exported {len(results_df)} results to CSV: {output_path}"
            )

        elif file_extension in ['.xlsx', '.xls']:
            results_df.to_excel(output_path, index=False, engine='openpyxl')
            logger.info(
                f"Exported {len(results_df)} results to Excel: {output_path}"
            )

        else:
            raise ValueError(
                f"Unsupported export format: '{file_extension}'. "
                "Supported formats: .csv, .xlsx, .xls"
            )

    except Exception as e:
        raise DataLoaderError(
            f"Failed to export results to '{output_path}': "
            f"{type(e).__name__} - {e}"
        ) from e

    # write error log
    if write_error_log and errors:
        error_log_path = output_path.with_suffix(
            output_path.suffix + '_errors.txt'
        )
        _write_error_log(errors, error_log_path, output_path)


def _write_error_log(
    errors: List[ParsingError],
    log_path: Path,
    source_path: Path
) -> None:
    """
    Error log file.

    Args:
        errors: Parsing errors
        log_path: Path of error log file
        source_path: Path of results file (for reference)
    """
    try:
        with open(log_path, 'w', encoding='utf-8') as f:
            # Header
            f.write("HLA Processing Error Log\n")
            f.write("=" * 50 + "\n")
            f.write(
                f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            )
            f.write(f"Results file: {source_path.name}\n")
            f.write(f"Total errors: {len(errors)}\n")
            f.write("=" * 50 + "\n\n")

            # group by type
            error_types = {}
            for error in errors:
                error_types[error.error_type] = (
                    error_types.get(error.error_type, 0) + 1
                )

            # summary
            f.write("Error Summary:\n")
            for error_type, count in sorted(error_types.items()):
                f.write(f"  - {error_type}: {count}\n")
            f.write("\n" + "-" * 50 + "\n\n")

            # details
            f.write("Detailed Errors:\n\n")
            for i, error in enumerate(errors, 1):
                f.write(f"Error {i}:\n")
                f.write(f"  Row ID: {error.row_id}\n")
                f.write(f"  Type: {error.error_type}\n")
                f.write(f"  Message: {error.message}\n")

                # error specifics
                if error.details:
                    f.write("  Details:\n")
                    for key, value in error.details.items():
                        f.write(f"    - {key}: {value}\n")

                f.write("\n")

        logger.info(
            f"Error log with {len(errors)} errors accessible here: {log_path}"
        )

    except Exception as e:
        logger.error(
            f"Failed to write error log to '{log_path}': "
            f"{type(e).__name__} - {e}"
        )
