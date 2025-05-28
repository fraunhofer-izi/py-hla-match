import unittest
import tempfile
import pandas as pd
from pathlib import Path
import logging

from py_hla_match.export import export_results
from py_hla_match.exceptions import DataLoaderError, ParsingError

logger = logging.getLogger(__name__)


class TestExport(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        self.sample_df = pd.DataFrame({
            'patient_id': [1, 2, 3],
            'donor_id': [4, 5, 6],
            'A_match_level1': ['ARD_MATCH', 'ALLELE_MISMATCH', 'ARD_MATCH'],
            'A_match_level2': ['ARD_MATCH', 'ARD_MATCH', 'ALLELE_MISMATCH'],
            'B_match_level1': ['ARD_MATCH', 'ARD_MATCH', 'ARD_MATCH'],
            'B_match_level2': ['ARD_MATCH', 'ARD_MATCH', 'ARD_MATCH']
        })

    def tearDown(self):
        self.temp_dir.cleanup()

    def create_parsing_errors(self):
        return [
            ParsingError(
                row_id=1,
                error_type='InvalidHLAString',
                message="Invalid HLA string 'INVALID' in column 'HLA_A1'",
                details={'column': 'HLA_A1', 'value': 'INVALID'}
            ),
            ParsingError(
                row_id=2,
                error_type='IncompleteLocus',
                message="Only one allele found for locus B",
                details={'locus': 'B'}
            )
        ]

    def test_export_to_csv(self):
        output_file = self.temp_path / "results.csv"
        export_results(self.sample_df, output_file)
        self.assertTrue(output_file.exists())
        loaded_df = pd.read_csv(output_file)
        self.assertEqual(len(loaded_df), 3)
        self.assertEqual(list(loaded_df.columns), list(self.sample_df.columns))

    def test_export_to_xlsx(self):
        output_file = self.temp_path / "results.xlsx"
        export_results(self.sample_df, output_file)
        self.assertTrue(output_file.exists())
        loaded_df = pd.read_excel(output_file)
        self.assertEqual(len(loaded_df), 3)
        self.assertEqual(list(loaded_df.columns), list(self.sample_df.columns))

    def test_export_to_xls(self):
        output_file = self.temp_path / "results.xls"
        export_results(self.sample_df, output_file)
        self.assertTrue(output_file.exists())
        loaded_df = pd.read_excel(output_file)
        self.assertEqual(len(loaded_df), 3)

    def test_export_creates_directory(self):
        output_file = self.temp_path / "subdir" / "nested" / "results.csv"
        export_results(self.sample_df, output_file)
        self.assertTrue(output_file.exists())
        self.assertTrue(output_file.parent.exists())

    def test_export_empty_dataframe(self):
        empty_df = pd.DataFrame()
        output_file = self.temp_path / "empty.csv"
        with self.assertLogs('py_hla_match.export', level='WARNING') as log:
            export_results(empty_df, output_file)
        self.assertIn("Exporting empty DataFrame", log.output[0])
        self.assertTrue(output_file.exists())

    def test_export_with_string_path(self):
        output_file = str(self.temp_path / "string_path.csv")
        export_results(self.sample_df, output_file)
        self.assertTrue(Path(output_file).exists())

    def test_export_unsupported_format(self):
        output_file = self.temp_path / "results.txt"
        with self.assertRaises(DataLoaderError) as context:
            export_results(self.sample_df, output_file)
        self.assertIn(
            "Unsupported export format: '.txt'", str(context.exception)
        )
        self.assertIn(
            "Supported formats: .csv, .xlsx, .xls", str(context.exception)
        )

    def test_export_with_error_log(self):
        output_file = self.temp_path / "results.csv"
        errors = self.create_parsing_errors()
        export_results(self.sample_df, output_file, errors=errors)
        # main file
        self.assertTrue(output_file.exists())
        # error log
        error_log = self.temp_path / "results.csv_errors.txt"
        self.assertTrue(error_log.exists())

    def test_error_log_content(self):
        output_file = self.temp_path / "results.xlsx"
        errors = self.create_parsing_errors()
        export_results(self.sample_df, output_file, errors=errors)
        error_log = self.temp_path / "results.xlsx_errors.txt"
        content = error_log.read_text()
        # header
        self.assertIn("HLA Processing Error Log", content)
        self.assertIn("Total errors: 2", content)
        # summary
        self.assertIn("InvalidHLAString: 1", content)
        self.assertIn("IncompleteLocus: 1", content)
        # details
        self.assertIn("Row ID: 1", content)
        self.assertIn("Row ID: 2", content)
        self.assertIn("Type: InvalidHLAString", content)
        self.assertIn("Type: IncompleteLocus", content)

    def test_export_without_error_log(self):
        output_file = self.temp_path / "results.csv"
        errors = self.create_parsing_errors()
        export_results(
            self.sample_df, output_file, errors=errors, write_error_log=False
        )
        # main file
        self.assertTrue(output_file.exists())
        # no error log
        error_log = self.temp_path / "results.csv_errors.txt"
        self.assertFalse(error_log.exists())

    def test_export_without_errors_no_log(self):
        output_file = self.temp_path / "results.csv"
        export_results(self.sample_df, output_file, errors=None)
        # main file
        self.assertTrue(output_file.exists())
        # no error(s) log
        error_log = self.temp_path / "results.csv_errors.txt"
        self.assertFalse(error_log.exists())

    def test_export_empty_errors_list_no_log(self):
        output_file = self.temp_path / "results.csv"
        export_results(self.sample_df, output_file, errors=[])
        # main file
        self.assertTrue(output_file.exists())
        # empty error log
        error_log = self.temp_path / "results.csv_errors.txt"
        self.assertFalse(error_log.exists())

    def test_error_log_with_error_details(self):
        output_file = self.temp_path / "results.csv"
        error_with_details = ParsingError(
            row_id=5,
            error_type='ComplexError',
            message='Complex error with multiple details',
            details={
                'column': 'HLA_A1',
                'value': 'BAD_VALUE',
                'additional_info': 'Extra context'
            }
        )
        export_results(
            self.sample_df, output_file, errors=[error_with_details]
        )
        error_log = self.temp_path / "results.csv_errors.txt"
        content = error_log.read_text()
        self.assertIn("column: HLA_A1", content)
        self.assertIn("value: BAD_VALUE", content)
        self.assertIn("additional_info: Extra context", content)

    def test_error_log_without_details(self):
        output_file = self.temp_path / "results.csv"
        error_without_details = ParsingError(
            row_id=10,
            error_type='SimpleError',
            message='Simple error without details',
            details=None
        )
        export_results(
            self.sample_df, output_file, errors=[error_without_details]
        )
        error_log = self.temp_path / "results.csv_errors.txt"
        content = error_log.read_text()
        self.assertIn("Row ID: 10", content)
        self.assertIn("Type: SimpleError", content)
        self.assertIn("Message: Simple error without details", content)
        # no details section
        self.assertNotIn(
            "Details:", content.split("Row ID: 10")[1].split("\n\n")[0]
        )

    def test_export_permission_error(self):
        # create read only file
        output_file = self.temp_path / "readonly.csv"
        output_file.touch()
        output_file.chmod(0o444)  # read-only
        try:
            with self.assertRaises(DataLoaderError) as context:
                export_results(self.sample_df, output_file)
            self.assertIn("Failed to export results", str(context.exception))
        except AssertionError:
            logger.warning(
                "Permission test may not work on this platform, skipping."
            )
            self.skipTest("Permission test may not work on this platform")
        finally:
            # restore permissions to delete file
            output_file.chmod(0o644)

    def test_export_logs_successful_export(self):
        output_file = self.temp_path / "results.csv"
        with self.assertLogs('py_hla_match.export', level='INFO') as log:
            export_results(self.sample_df, output_file)
        self.assertTrue(
            any("Exported 3 results to CSV" in msg for msg in log.output)
        )

    def test_export_error_log_creation_logged(self):
        output_file = self.temp_path / "results.csv"
        errors = self.create_parsing_errors()
        with self.assertLogs('py_hla_match.export', level='INFO') as log:
            export_results(self.sample_df, output_file, errors=errors)
        self.assertTrue(
            any(
                "Error log with 2 errors accessible here" in msg for msg
                in log.output
            )
        )


if __name__ == "__main__":
    unittest.main()
