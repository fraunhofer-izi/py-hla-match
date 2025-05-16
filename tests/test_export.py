import os
import unittest
import tempfile

import pandas as pd

from py_hla_match.export import PairwiseMatchResult, BestMatchResult
from py_hla_match.parser import HLADataSource


class TestExport(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up file paths for test cases."""
        TEST_DIR_PATH = os.path.dirname(os.path.realpath(__file__))
        cls.valid_csv = os.path.join(
            TEST_DIR_PATH, "resources", "hla_test_data.csv"
        )
        cls.valid_excel = os.path.join(
            TEST_DIR_PATH, "resources", "hla_test_data.xlsx"
        )
        # temporary output directory
        cls.temp_dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        """Clean up temporary directory."""
        # Remove temp directory if it exists
        import shutil
        if os.path.exists(cls.temp_dir):
            shutil.rmtree(cls.temp_dir)

    def tearDown(self):
        """Clean up temporary files created during tests."""
        for file in ["match_results.csv", "best_match_result.csv"]:
            if os.path.exists(file):
                os.remove(file)

    def test_pairwise_match_object_based_basic_resolution(self):
        """
        Test pairwise matching with object-based data and basic resolution.
        """
        # basic parsing to get a list of individuals
        individuals = HLADataSource(self.valid_csv).parse()
        # create pairwise match object
        match_result = PairwiseMatchResult(
            individuals, individuals, resolution="basic"
        ).to_df()
        # check that the result is a pandas DataFrame
        self.assertIsInstance(match_result, pd.DataFrame)
        # check that the result has the expected number of rows
        self.assertEqual(len(match_result), 8)
        # check that the result has the expected number of columns
        self.assertEqual(len(match_result.columns), 5)
        # check that the resulting dataframe has correct column names
        expected_columns = ["A", "B", "C", "DPA1", "DRB1"]
        self.assertEqual(list(match_result.columns), expected_columns)

    def test_best_match_object_based_basic_resolution(self):
        """Test best match with object-based data and basic resolution."""
        # basic parsing to get a list of individuals
        individuals = HLADataSource(self.valid_csv).parse()
        # create pairwise match object
        match_result = BestMatchResult(
            individuals[3], individuals, resolution="basic"
        )
        best_match = match_result.get_best_match()
        best_match_idx = match_result.get_best_match_target_idx()
        # third individual in the list should be the best match
        self.assertEqual(best_match_idx, 3)
        self.assertEqual(best_match, individuals[3])

    def test_pairwise_match_file_based(self):
        """Test pairwise matching with file-based data."""
        individuals = HLADataSource(self.valid_csv).parse()
        match_result = PairwiseMatchResult(individuals, individuals)
        df = match_result.to_df()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 8)  # !! with 8 individuals in test file

    def test_pairwise_match_excel_input(self):
        """Test pairwise matching with Excel input file."""
        individuals = HLADataSource(self.valid_csv).parse()
        match_result = PairwiseMatchResult(individuals, individuals)
        df = match_result.to_df()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df.columns), 5)  # !! with 5 loci in test file

    def test_pairwise_match_high_resolution(self):
        """Test pairwise matching with high resolution."""
        individuals = HLADataSource(self.valid_csv).parse()
        match_result = PairwiseMatchResult(
            individuals, individuals, resolution="high"
        ).to_df()
        self.assertIsInstance(match_result, pd.DataFrame)
        # !! with 5 loci intest file
        self.assertEqual(len(match_result.columns), 5)

    def test_pairwise_match_full_resolution(self):
        """Not implemented yet."""
        pass

    def test_export_to_csv(self):
        """Test exporting results to CSV."""
        individuals = HLADataSource(self.valid_csv).parse()
        match_result = PairwiseMatchResult(individuals, individuals)
        output_file = os.path.join(self.temp_dir, "test_output.csv")
        match_result.to_csv(output_file)
        self.assertTrue(os.path.exists(output_file))
        df = pd.read_csv(output_file)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 8)  # !! with 8 individuals in test file

    def test_export_to_excel(self):
        """Test exporting results to Excel."""
        individuals = HLADataSource(self.valid_csv).parse()
        match_result = PairwiseMatchResult(individuals, individuals)
        output_file = os.path.join(self.temp_dir, "test_output.xlsx")
        match_result.to_excel(output_file)
        self.assertTrue(os.path.exists(output_file))
        # Verify file content
        df = pd.read_excel(output_file)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 8)  # !! with 8 individuals in test file

    def test_invalid_resolution(self):
        """Test that invalid resolution raises ValueError."""
        individuals = HLADataSource(self.valid_csv).parse()
        with self.assertRaises(ValueError):
            PairwiseMatchResult(individuals, individuals, resolution="invalid")

    def test_unsupported_file_format(self):
        """Test that unsupported file format raises ValueError."""
        unsupported_file = os.path.join(self.temp_dir, "test.txt")
        with open(unsupported_file, 'w') as f:
            f.write("test")
        with self.assertRaises(ValueError):
            PairwiseMatchResult(unsupported_file, self.valid_csv)
