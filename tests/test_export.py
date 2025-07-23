import os
import unittest
import tempfile

import pandas as pd

from py_hla_match.export import PairwiseMatch
from py_hla_match.parser import HLADataSource


def _generate_tmp_file_source_and_target(source_df: pd.DataFrame,
                                         target_df: pd.DataFrame) -> tuple[HLADataSource, HLADataSource]:
    """
    Generate temporary files for source and target HLA data.
    """
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as sf:
        source_df.to_csv(sf.name, index=False)
        source_path = sf.name
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tf:
        target_df.to_csv(tf.name, index=False)
        target_path = tf.name

    source = HLADataSource(source_path)
    target = HLADataSource(target_path)

    return source, target


class TestExport(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up file paths for test cases."""
        TEST_DIR_PATH = os.path.dirname(os.path.realpath(__file__))
        cls.valid_csv = os.path.join(
            TEST_DIR_PATH, "resources", "hla_test_data.csv"
        )
        cls.invalid_csv = os.path.join(
            TEST_DIR_PATH, "resources", "hla_test_data_malformed.csv"
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
        for file in ["match_results.csv"]:
            if os.path.exists(file):
                os.remove(file)

    def test_valid_csv_no_streaming(self):
        """Test parsing a valid CSV file without streaming."""
        source = HLADataSource(self.valid_csv)
        target = HLADataSource(self.valid_csv)
        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(self.temp_dir, "match_results.csv"),
            resolution="basic",
            stream=False
        )
        # assert that the generated file exists
        self.assertTrue(os.path.exists(pairwise_match.result_file))
        # assert content of DataFrame
        result_df = pairwise_match.to_df()
        self.assertIsInstance(result_df, pd.DataFrame)
        self.assertEqual(result_df.shape[0], 8)
        self.assertEqual(result_df.shape[1], 5)

    def test_valid_csv_streaming(self):
        """Test parsing a valid CSV file with streaming."""
        source = HLADataSource(self.valid_csv)
        target = HLADataSource(self.valid_csv)
        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(self.temp_dir, "match_results.csv"),
            resolution="basic",
            stream=True,
            chunk_size=2
        )
        # assert that the generated file exists
        self.assertTrue(os.path.exists(pairwise_match.result_file))
        # assert exception gets thrown when trying to access the streamed df
        with self.assertRaises(ValueError):
            pairwise_match.to_df()

    def test_valid_excel_streaming(self):
        """Test parsing a valid Excel file with streaming."""
        source = HLADataSource(self.valid_excel)
        target = HLADataSource(self.valid_excel)
        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(self.temp_dir, "match_results_excel_streaming.csv"),
            resolution="basic",
            stream=True,
            chunk_size=2
        )
        # Assert that the generated file exists
        self.assertTrue(os.path.exists(pairwise_match.result_file))
        # Assert exception gets thrown when trying to access the streamed DataFrame
        with self.assertRaises(ValueError):
            pairwise_match.to_df()

    def test_valid_excel_no_streaming(self):
        """Test parsing a valid Excel file without streaming."""
        source = HLADataSource(self.valid_excel)
        target = HLADataSource(self.valid_excel)
        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(self.temp_dir, "match_results_excel_no_streaming.csv"),
            resolution="basic",
            stream=False,
            chunk_size=2
        )
        # Assert that the result is stored in memory as a DataFrame
        self.assertIsInstance(pairwise_match.result, pd.DataFrame)
        self.assertGreater(len(pairwise_match.result), 0)  # Ensure rows are present
        self.assertGreater(len(pairwise_match.result.columns), 0)  # Ensure loci columns are present

    def test_invalid_csv_streaming(self):
        """Test parsing an invalid CSV file with streaming."""
        source = HLADataSource(self.invalid_csv)
        target = HLADataSource(self.invalid_csv)
        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(self.temp_dir, "match_results_invalid_csv_streaming.csv"),
            resolution="basic",
            stream=True,
            chunk_size=2
        )
        # Assert that the generated file exists
        self.assertTrue(os.path.exists(pairwise_match.result_file))
        # Assert that the file contains valid rows despite malformed data
        df = pd.read_csv(pairwise_match.result_file)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)  # Ensure rows are written
        self.assertGreater(len(df.columns), 0)  # Ensure loci columns are present

    def test_invalid_csv_no_streaming(self):
        """Test parsing an invalid CSV file without streaming."""
        source = HLADataSource(self.invalid_csv)
        target = HLADataSource(self.invalid_csv)
        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(self.temp_dir, "match_results_invalid_csv_no_streaming.csv"),
            resolution="basic",
            stream=False,
            chunk_size=2
        )
        # Assert that the result is stored in memory as a DataFrame
        self.assertIsInstance(pairwise_match.result, pd.DataFrame)
        self.assertGreater(len(pairwise_match.result), 0)  # Ensure rows are present
        self.assertGreater(len(pairwise_match.result.columns), 0)  # Ensure loci columns are present
        # Assert that errors were logged for malformed data
        with self.assertLogs("py_hla_match.parser", level="ERROR") as log_context:
            source.parse(stream=False)
            error_logs = [record for record in log_context.output if "Encountered malformed HLA String" in record]
            self.assertGreater(len(error_logs), 0, "Expected malformed HLA string log entries")

    def test_length_mismatch_raises(self):
        """
        PairwiseMatchResult must raise a ValueError.
        """
        source_df = pd.DataFrame(
            {"A1": ["A*01:01", "A*01:01"],
             "A2": ["A*02:01", "A*02:01"]}
        )  # 2 rows
        target_df = pd.DataFrame(
            {"A1": ["A*01:01", "A*01:01", "A*01:01"],
             "A2": ["A*02:01", "A*02:01", "A*02:01"]}
        )  # 3 rows

        source, target = _generate_tmp_file_source_and_target(source_df, target_df)

        with self.assertRaises(ValueError):
            PairwiseMatch(
                source=source,
                target=target,
                storage_filename=os.path.join(
                    self.__class__.temp_dir, "len_mismatch.csv"
                ),
                resolution="basic",
                stream=False
            )

    def test_unexpected_locus_raises(self):
        """
        Patient introduces a third locus in second row (by mistake). Export
        raises ValueError.
        """
        source_df = pd.DataFrame(
            {
                "A1": ["A*01:01", "A*01:01"],
                "A2": ["A*02:01", "A*02:01"],
                "B1": ["B*07:183", "C*03:02:06"],
                "B2": ["B*07:183", "C*03:02:06"]
             }
        )
        target_df = pd.DataFrame(
            {
                "A1": ["A*01:01", "A*01:01"],
                "A2": ["A*02:01", "A*02:01"],
                "B1": ["B*07:183", "B*07:183"],
                "B2": ["B*07:183", "B*07:183"]
             }
        )

        source, target = _generate_tmp_file_source_and_target(source_df, target_df)

        with self.assertRaises(ValueError):
            PairwiseMatch(
                source=source,
                target=target,
                storage_filename=os.path.join(
                    self.__class__.temp_dir, "bad_locus.csv"
                ),
                resolution="basic",
                stream=False
            )

