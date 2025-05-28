import unittest
import pandas as pd
from py_hla_match.parser import HLAParser
from py_hla_match.models import Individual


class TestHLAParser(unittest.TestCase):

    def setUp(self):
        self.parser = HLAParser()

    def create_single_df(self):
        return pd.DataFrame({
            'patient_id': [1, 2, 3],
            'HLA_A1': ['A*02:01', 'A*24:02', 'A*11:01'],
            'HLA_A2': ['A*03:01', 'A*25:01', 'A*31:01'],
            'HLA_B1': ['B*07:02', 'B*13:02', 'B*35:01'],
            'HLA_B2': ['B*08:01', 'B*44:03', 'B*51:01']
        })

    def create_paired_df(self):
        return pd.DataFrame({
            'pair_id': [1, 2],
            'pat_A1': ['A*02:01', 'A*24:02'],
            'pat_A2': ['A*03:01', 'A*25:01'],
            'pat_B1': ['B*07:02', 'B*13:02'],
            'pat_B2': ['B*08:01', 'B*44:03'],
            'don_A1': ['A*02:01', 'A*11:01'],
            'don_A2': ['A*03:01', 'A*31:01'],
            'don_B1': ['B*07:02', 'B*35:01'],
            'don_B2': ['B*08:01', 'B*51:01']
        })

    def create_panel_df(self):
        return pd.DataFrame({
            'family_id': [1],
            'mother_A1': ['A*02:01'],
            'mother_A2': ['A*03:01'],
            'father_A1': ['A*24:02'],
            'father_A2': ['A*25:01'],
            'child_A1': ['A*02:01'],
            'child_A2': ['A*24:02']
        })

    def test_parse_single_structure(self):
        df = self.create_single_df()
        column_mapping = {
            'individual': ['HLA_A1', 'HLA_A2', 'HLA_B1', 'HLA_B2']
        }
        results = self.parser.parse(
            df, 'single', column_mapping, id_column='patient_id'
        )
        self.assertEqual(len(results), 3)
        self.assertEqual(results[0][0], 1)  # ID
        self.assertIsInstance(results[0][1], Individual)

    def test_parse_paired_structure(self):
        df = self.create_paired_df()
        column_mapping = {
            'patient': ['pat_A1', 'pat_A2', 'pat_B1', 'pat_B2'],
            'donor': ['don_A1', 'don_A2', 'don_B1', 'don_B2']
        }
        results = self.parser.parse(
            df, 'paired', column_mapping, id_column='pair_id'
        )
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0][0], 1)  # ID
        self.assertIsInstance(results[0][1], Individual)  # Patient
        self.assertIsInstance(results[0][2], Individual)  # Donor

    def test_parse_panel_structure(self):
        df = self.create_panel_df()
        column_mapping = {
            'mother': ['mother_A1', 'mother_A2'],
            'father': ['father_A1', 'father_A2'],
            'child': ['child_A1', 'child_A2']
        }
        results = self.parser.parse(
            df, 'panel', column_mapping, id_column='family_id'
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0][0], 1)  # ID
        self.assertIsInstance(results[0][1], dict)
        self.assertIn('mother', results[0][1])
        self.assertIn('father', results[0][1])
        self.assertIn('child', results[0][1])

    def test_parse_without_id_column(self):
        df = pd.DataFrame({
            'HLA_A1': ['A*02:01', 'A*24:02'],
            'HLA_A2': ['A*03:01', 'A*25:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        results = self.parser.parse(df, 'single', column_mapping)
        self.assertEqual(results[0][0], 0)  # First row index
        self.assertEqual(results[1][0], 1)  # Second row index

    def test_parse_with_row_index_column(self):
        df = pd.DataFrame({
            '_row_index': [10, 20],
            'HLA_A1': ['A*02:01', 'A*24:02'],
            'HLA_A2': ['A*03:01', 'A*25:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        results = self.parser.parse(df, 'single', column_mapping)
        self.assertEqual(results[0][0], 10)
        self.assertEqual(results[1][0], 20)

    def test_parse_invalid_hla_string_creates_error(self):
        df = pd.DataFrame({
            'id': [1],
            'HLA_A1': ['INVALID*STRING'],
            'HLA_A2': ['A*03:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        results = self.parser.parse(df, 'single', column_mapping, 'id')
        self.assertEqual(len(results), 1)
        # creates 2 errors: InvalidHLAString + IncompleteLocus
        self.assertEqual(len(self.parser.errors), 2)
        error_types = [e.error_type for e in self.parser.errors]
        self.assertIn('InvalidHLAString', error_types)
        self.assertIn('IncompleteLocus', error_types)

    def test_parse_incomplete_locus_creates_error(self):
        df = pd.DataFrame({
            'id': [1],
            'HLA_A1': ['A*02:01'],
            'HLA_B1': ['B*07:02']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_B1']}
        self.parser.parse(df, 'single', column_mapping, 'id')
        self.assertEqual(len(self.parser.errors), 2)
        for error in self.parser.errors:
            self.assertEqual(error.error_type, 'IncompleteLocus')

    def test_parse_too_many_alleles_creates_error(self):
        df = pd.DataFrame({
            'id': [1],
            'HLA_A1': ['A*02:01'],
            'HLA_A2': ['A*03:01'],
            'HLA_A3': ['A*11:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2', 'HLA_A3']}
        self.parser.parse(df, 'single', column_mapping, 'id')
        # catch exception and log as ParseError
        self.assertEqual(len(self.parser.errors), 1)
        self.assertEqual(self.parser.errors[0].error_type, 'ParseError')
        self.assertIn(
            "Found 3 alleles for locus A", self.parser.errors[0].message
        )

    def test_parse_empty_values_are_skipped(self):
        df = pd.DataFrame({
            'id': [1],
            'HLA_A1': ['A*02:01'],
            'HLA_A2': [''],  # Empty string
            'HLA_B1': ['B*07:02'],
            'HLA_B2': ['B*08:01']
        })
        column_mapping = {
            'individual': ['HLA_A1', 'HLA_A2', 'HLA_B1', 'HLA_B2']
        }
        results = self.parser.parse(df, 'single', column_mapping, 'id')
        individual = results[0][1]
        self.assertEqual(len(individual.hla_data), 2)  # A and B loci

    def test_invalid_structure_raises_error(self):
        df = self.create_single_df()
        with self.assertRaises(ValueError) as context:
            self.parser.parse(df, 'invalid', {}, 'patient_id')
        self.assertIn("Invalid structure 'invalid'", str(context.exception))

    def test_missing_required_keys_for_single(self):
        df = self.create_single_df()
        with self.assertRaises(ValueError) as context:
            self.parser.parse(df, 'single', {'wrong_key': []}, 'patient_id')
        self.assertIn("column_mapping must contain 'individual' key",
                      str(context.exception))

    def test_missing_required_keys_for_paired(self):
        df = self.create_paired_df()
        with self.assertRaises(ValueError) as context:
            self.parser.parse(df, 'paired', {'patient': []}, 'pair_id')
        self.assertIn("column_mapping must contain 'patient' and 'donor' keys",
                      str(context.exception))

    def test_nonexistent_columns_in_mapping(self):
        df = self.create_single_df()
        column_mapping = {'individual': ['MISSING_COL1', 'MISSING_COL2']}
        with self.assertRaises(ValueError) as context:
            self.parser.parse(df, 'single', column_mapping, 'patient_id')
        self.assertIn("Columns not found in DataFrame", str(context.exception))

    def test_nonexistent_id_column(self):
        df = self.create_single_df()
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        with self.assertRaises(ValueError) as context:
            self.parser.parse(df, 'single', column_mapping, 'missing_id')
        self.assertIn(
            "ID column 'missing_id' not found", str(context.exception)
        )

    def test_get_error_summary_no_errors(self):
        summary = self.parser.get_error_summary()
        self.assertEqual(summary, "No parsing errors encountered.")

    def test_get_error_summary_with_errors(self):
        df = pd.DataFrame({
            'id': [1, 2],
            'HLA_A1': ['INVALID1', 'INVALID2'],
            'HLA_A2': ['A*03:01', 'A*25:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        self.parser.parse(df, 'single', column_mapping, 'id')
        summary = self.parser.get_error_summary()
        # 2 errors per invalid string
        self.assertIn("Total errors: 4", summary)
        self.assertIn("InvalidHLAString: 2", summary)
        self.assertIn("IncompleteLocus: 2", summary)

    def test_errors_reset_between_parses(self):
        df1 = pd.DataFrame({
            'id': [1],
            'HLA_A1': ['INVALID'],
            'HLA_A2': ['A*03:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        self.parser.parse(df1, 'single', column_mapping, 'id')
        # 2 errors per invalid string
        self.assertEqual(len(self.parser.errors), 2)

        df2 = pd.DataFrame({
            'id': [1],
            'HLA_A1': ['A*02:01'],
            'HLA_A2': ['A*03:01']
        })
        self.parser.parse(df2, 'single', column_mapping, 'id')
        self.assertEqual(len(self.parser.errors), 0)

    def test_parse_continues_after_row_error(self):
        df = pd.DataFrame({
            'id': [1, 2, 3],
            'HLA_A1': ['A*02:01', 'TOTALLY_INVALID', 'A*11:01'],
            'HLA_A2': ['A*03:01', 'ALSO_INVALID', 'A*31:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        results = self.parser.parse(df, 'single', column_mapping, 'id')
        # failed row should not stop parsing
        self.assertEqual(len(results), 3)
        # first and third allele valid
        self.assertEqual(len(results[0][1].hla_data), 1)
        self.assertEqual(len(results[2][1].hla_data), 1)
        # no valid data in second row
        self.assertEqual(len(results[1][1].hla_data), 0)

    def test_error_details_are_recorded(self):
        df = pd.DataFrame({
            'id': [1],
            'HLA_A1': ['INVALID*STRING'],
            'HLA_A2': ['A*03:01']
        })
        column_mapping = {'individual': ['HLA_A1', 'HLA_A2']}
        self.parser.parse(df, 'single', column_mapping, 'id')
        # InvalidHLAString
        invalid_error = [e for e in self.parser.errors
                         if e.error_type == 'InvalidHLAString'][0]
        self.assertEqual(invalid_error.row_id, 1)
        self.assertIn('HLA_A1', invalid_error.details['column'])
        self.assertEqual(invalid_error.details['value'], 'INVALID*STRING')


if __name__ == "__main__":
    unittest.main()
