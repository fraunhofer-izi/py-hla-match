import logging  # NOQA
import unittest
import os
import tempfile
import pandas as pd
from py_hla_match.models import Individual
from py_hla_match.parser import HLADataSource


class TestParser(unittest.TestCase):

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
        cls.invalid_csv = os.path.join(
            TEST_DIR_PATH, "resources", "hla_test_data_malformed.csv"
        )

    def _tmp_csv(self, df: pd.DataFrame) -> str:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as f:
            df.to_csv(f.name, index=False)
            return f.name

    def test_parse_valid_csv(self):
        """Test parsing a valid CSV file."""
        parser = HLADataSource(self.valid_csv)
        individuals = parser.parse()
        self.assertIsInstance(individuals, list)
        self.assertEqual(len(individuals), 8)
        self.assertTrue(
            all(isinstance(ind, Individual) for ind in individuals)
        )
        self.assertEqual("DRB1", individuals[0].hla_data[4].locus)

    def test_parse_malformed_hla_strings_get_logged(self):
        """
        Test parsing an invalid CSV file logs MalformedHLAStringError twice.
        """
        parser = HLADataSource(self.invalid_csv)
        with self.assertLogs(
                "py_hla_match.parser",
                level="ERROR"
        ) as log_context:
            parser.parse()
        error_logs = [
            record for record in log_context.output if
            "Encountered malformed HLA String" in record
        ]
        self.assertEqual(
            len(error_logs),
            3,
            "Expected exactly two malformed HLA string log entries"
        )

    def test_parse_valid_excel(self):
        """Test parsing a valid Excel file."""
        parser = HLADataSource(self.valid_excel)
        individuals = parser.parse()
        self.assertIsInstance(individuals, list)
        self.assertEqual(len(individuals), 8)
        self.assertTrue(
            all(isinstance(ind, Individual) for ind in individuals)
        )
        self.assertEqual("DRB1", individuals[0].hla_data[4].locus)

    def test_stream_csv(self):
        """Test streaming parsing of a valid CSV file."""
        parser = HLADataSource(self.valid_csv)
        individuals = parser.parse(stream=True, chunk_size=2)
        n_streamed = 0
        for individual in individuals:
            n_streamed += 1
            self.assertIsInstance(individual, Individual)
        self.assertEqual(n_streamed, 8)

    def test_stream_csv_with_col_indices(self):
        """Test streaming parsing of a valid CSV file with column indices."""
        parser = HLADataSource(self.valid_csv, col_idx_start=2, col_idx_stop=5)
        individuals = parser.parse(stream=True, chunk_size=2)
        n_streamed = 0
        for individual in individuals:
            n_streamed += 1
            self.assertIsInstance(individual, Individual)
            self.assertEqual(len(individual.hla_data), 2)
        self.assertEqual(n_streamed, 8)

    def test_stream_excel(self):
        """Test streaming parsing of a valid Excel file."""
        parser = HLADataSource(self.valid_excel)
        individuals = parser.parse(stream=True, chunk_size=2)
        n_streamed = 0
        for individual in individuals:
            n_streamed += 1
            self.assertIsInstance(individual, Individual)
        self.assertEqual(n_streamed, 8)

    def test_stream_excel_with_col_indices(self):
        """Test streaming parsing of a valid Excel file with column indices."""
        parser = HLADataSource(
            self.valid_excel, col_idx_start=2, col_idx_stop=5
            )
        individuals = parser.parse(stream=True, chunk_size=2)
        n_streamed = 0
        for individual in individuals:
            n_streamed += 1
            self.assertIsInstance(individual, Individual)
            self.assertEqual(len(individual.hla_data), 2)
        self.assertEqual(n_streamed, 8)

    def test_non_stream_parsing_handles_zero_start_index(self):
        """
        Test non-stream correctly parses slices when col_idx_start is 0.
        """
        parser = HLADataSource(self.valid_csv, col_idx_start=0, col_idx_stop=3)

        individuals = parser.parse()

        self.assertEqual(
            len(individuals[0].hla_data),
            2,
            "Parser failed to slice correctly when start index was 0"
        )

    # NOTE: this test failes due to the same bug (col_idx_start is 0)
    def test_slice_consistency_non_stream_vs_stream(self):
        """
        Tests that slicing is consistent between stream and non-stream modes.
        """
        df = pd.DataFrame(
            [["A*01:01", "A*02:01", "B*07:02", "B*08:01"]],
            columns=["A1", "A2", "B1", "B2"],
        )
        path = self._tmp_csv(df)

        non_stream_parser = HLADataSource(
            path, col_idx_start=0, col_idx_stop=1
        )
        stream_parser = HLADataSource(
            path, col_idx_start=0, col_idx_stop=1
        )

        non_stream_individuals = non_stream_parser.parse()
        stream_individuals = list(
            stream_parser.parse(stream=True, chunk_size=1)
        )
        self.assertEqual(
            len(non_stream_individuals[0].hla_data),
            len(stream_individuals[0].hla_data),
            "Non-streaming and streaming result has different number of pairs."
        )


if __name__ == "__main__":
    unittest.main()
