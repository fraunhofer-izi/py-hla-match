import unittest
import tempfile
import pandas as pd
from pathlib import Path

from py_hla_match.loader import HLADataLoader
from py_hla_match.exceptions import (
    FileNotFoundError,
    DataLoaderError,
    EmptyDataError,
    UnsupportedFileTypeError
)


class TestHLADataLoader(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def create_test_csv(self, filename, content):
        filepath = self.temp_path / filename
        with open(filepath, 'w') as f:
            f.write(content)
        return filepath

    def create_test_excel(self, filename, data_dict):
        filepath = self.temp_path / filename
        df = pd.DataFrame(data_dict)
        df.to_excel(filepath, index=False)
        return filepath

    def test_load_valid_csv(self):
        content = \
            "donor_id,HLA_A1,HLA_A2\n1,A*02:01,A*03:01\n2,A*24:02,A*25:01"
        filepath = self.create_test_csv("test.csv", content)
        df = HLADataLoader.load(filepath)
        self.assertEqual(len(df), 2)
        self.assertEqual(
            list(df.columns), ['donor_id', 'HLA_A1', 'HLA_A2', '_row_index']
        )

    def test_load_csv_with_specific_columns(self):
        content = \
            "donor_id,HLA_A1,HLA_A2,HLA_B1,HLA_B2\n1,A*02:01,A*03:01," \
            "B*07:02,B*08:01"
        filepath = self.create_test_csv("test.csv", content)
        df = HLADataLoader.load(
            filepath, columns=['donor_id', 'HLA_A1', 'HLA_A2']
        )
        self.assertEqual(len(df.columns), 4)  # includes _row_index
        self.assertIn('HLA_A1', df.columns)
        self.assertNotIn('HLA_B1', df.columns)

    def test_load_csv_with_id_column(self):
        content = \
            "donor_id,HLA_A1,HLA_A2\n1,A*02:01,A*03:01\n2,A*24:02,A*25:01"
        filepath = self.create_test_csv("test.csv", content)
        df = HLADataLoader.load(filepath, id_column='donor_id')
        self.assertNotIn('_row_index', df.columns)
        self.assertIn('donor_id', df.columns)

    def test_load_csv_with_custom_delimiter(self):
        content = \
            "donor_id;HLA_A1;HLA_A2\n1;A*02:01;A*03:01\n2;A*24:02;A*25:01"
        filepath = self.create_test_csv("test.csv", content)
        df = HLADataLoader.load(filepath, delimiter=';')
        self.assertEqual(len(df), 2)
        self.assertEqual(df.iloc[0]['HLA_A1'], 'A*02:01')

    def test_load_excel_xlsx(self):
        data = {
            'patient_id': [1, 2],
            'HLA_A1': ['A*02:01', 'A*24:02'],
            'HLA_A2': ['A*03:01', 'A*25:01']
        }
        filepath = self.create_test_excel("test.xlsx", data)
        df = HLADataLoader.load(filepath)
        self.assertEqual(len(df), 2)
        self.assertEqual(df.iloc[0]['HLA_A1'], 'A*02:01')

    def test_load_excel_with_sheet_name(self):
        filepath = self.temp_path / "test_sheets.xlsx"
        with pd.ExcelWriter(filepath) as writer:
            pd.DataFrame({'A': [1, 2]}).to_excel(
                writer, sheet_name='Sheet1', index=False
            )
            pd.DataFrame({'B': [3, 4]}).to_excel(
                writer, sheet_name='Sheet2', index=False
            )
        df = HLADataLoader.load(filepath, sheet_name='Sheet2')
        self.assertIn('B', df.columns)
        self.assertNotIn('A', df.columns)

    def test_load_excel_with_sheet_index(self):
        filepath = self.temp_path / "test_sheets.xlsx"
        with pd.ExcelWriter(filepath) as writer:
            pd.DataFrame({'A': [1, 2]}).to_excel(
                writer, sheet_name='First', index=False
            )
            pd.DataFrame({'B': [3, 4]}).to_excel(
                writer, sheet_name='Second', index=False
            )
        df = HLADataLoader.load(filepath, sheet_name=1)
        self.assertIn('B', df.columns)
        self.assertNotIn('A', df.columns)

    def test_load_file_not_found(self):
        filepath = self.temp_path / "nonexistent.csv"
        with self.assertRaises(FileNotFoundError):
            HLADataLoader.load(filepath)

    def test_load_directory_instead_of_file(self):
        with self.assertRaises(FileNotFoundError):
            HLADataLoader.load(self.temp_path)

    def test_load_empty_file(self):
        filepath = self.create_test_csv("empty.csv", "")
        with self.assertRaises(EmptyDataError):
            HLADataLoader.load(filepath)

    def test_load_csv_with_only_headers(self):
        content = "donor_id,HLA_A1,HLA_A2"
        filepath = self.create_test_csv("headers_only.csv", content)
        with self.assertRaises(EmptyDataError):
            HLADataLoader.load(filepath)

    def test_load_unsupported_file_type(self):
        filepath = self.temp_path / "test.txt"
        filepath.write_text("some content")
        with self.assertRaises(UnsupportedFileTypeError):
            HLADataLoader.load(filepath)

    def test_load_with_nonexistent_id_column(self):
        content = "donor_id,HLA_A1,HLA_A2\n1,A*02:01,A*03:01"
        filepath = self.create_test_csv("test.csv", content)
        with self.assertRaises(DataLoaderError) as context:
            HLADataLoader.load(filepath, id_column='missing_column')
        self.assertIn(
            "Specified ID column 'missing_column' not found",
            str(context.exception)
        )

    def test_load_with_nonexistent_columns(self):
        content = "donor_id,HLA_A1,HLA_A2\n1,A*02:01,A*03:01"
        filepath = self.create_test_csv("test.csv", content)
        with self.assertRaises(DataLoaderError) as context:
            HLADataLoader.load(filepath, columns=['missing1', 'missing2'])
        self.assertIn("Usecols do not match columns", str(context.exception))

    def test_load_with_partially_matching_columns(self):
        content = "donor_id,HLA_A1,HLA_A2\n1,A*02:01,A*03:01"
        filepath = self.create_test_csv("test.csv", content)
        with self.assertRaises(DataLoaderError) as context:
            HLADataLoader.load(
                filepath, columns=['donor_id', 'missing_column']
            )
        self.assertIn("Usecols do not match columns", str(context.exception))

    def test_all_data_loaded_as_strings(self):
        content = "id,HLA_A1,count\n1,A*02:01,123\n2,A*24:02,456"
        filepath = self.create_test_csv("test.csv", content)
        df = HLADataLoader.load(filepath)
        self.assertEqual(df['count'].dtype, object)  # String type in pandas
        self.assertEqual(df.iloc[0]['count'], '123')  # Not int

    def test_no_nan_conversion(self):
        content = "id,HLA_A1,HLA_A2\n1,A*02:01,\n2,,A*24:02"
        filepath = self.create_test_csv("test.csv", content)
        df = HLADataLoader.load(filepath)
        self.assertEqual(df.iloc[0]['HLA_A2'], '')  # Empty string, not NaN
        self.assertEqual(df.iloc[1]['HLA_A1'], '')  # Empty string, not NaN

    def test_row_index_sequential(self):
        content = \
            "HLA_A1,HLA_A2\nA*02:01,A*03:01\nA*24:02,A*25:01\nA*11:01,A*31:01"
        filepath = self.create_test_csv("test.csv", content)
        df = HLADataLoader.load(filepath)
        self.assertEqual(list(df['_row_index']), [0, 1, 2])

    def test_validate_file_exists(self):
        filepath = self.temp_path / "nonexistent.csv"
        with self.assertRaises(FileNotFoundError):
            HLADataLoader._validate_file(filepath)

    def test_validate_file_is_directory(self):
        with self.assertRaises(FileNotFoundError):
            HLADataLoader._validate_file(self.temp_path)

    def test_validate_file_is_empty(self):
        filepath = self.temp_path / "empty.txt"
        filepath.touch()  # Creates empty file
        with self.assertRaises(EmptyDataError):
            HLADataLoader._validate_file(filepath)


if __name__ == "__main__":
    unittest.main()
