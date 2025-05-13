import os
import unittest

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

    def test_pairwise_match_object_based_basic_resolution(self):
        """Test pairwise matching with object-based data and basic resolution."""
        # basic parsing to get a list of individuals
        individuals = HLADataSource(self.valid_csv).parse()
        # create pairwise match object
        match_result = PairwiseMatchResult(individuals, individuals, resolution="basic").to_df()
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
        match_result = BestMatchResult(individuals[3], individuals, resolution="basic")
        best_match = match_result.get_best_match()
        best_match_idx = match_result.get_best_match_target_idx()
        # third individual in the list should be the best match
        self.assertEqual(best_match_idx, 3)
        self.assertEqual(best_match, individuals[3])

