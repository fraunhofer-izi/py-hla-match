import logging  # NOQA
import unittest
import os
import tempfile
import csv
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
        Test parsing an invalid CSV  file with three incorrect fields logs
        MalformedHLAStringError exactly three times.
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
            "Expected exactly three malformed HLA string log entries"
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

    def test_row_idx_start_respected_csv(self):
        """
        Test that row_idx_start is respected in CSV parsing (stream and
        non-stream).
        """
        # create a CSV with 2 header/garbage rows
        # row 0: metadata
        # row 1: columns
        # row 2: actual data
        data = [
            ["MetaData", "Version 1.0"],
            ["A_1", "A_2"],
            ["A*01:01", "A*02:01"]
        ]

        # write raw CSV
        with tempfile.NamedTemporaryFile(
            delete=False, suffix=".csv", mode='w', newline=''
        ) as f:
            writer = csv.writer(f)
            writer.writerows(data)
            tmp_path = f.name

        try:
            # configure parser to start at row 2 (skipping 0 and 1)
            parser = HLADataSource(tmp_path, row_idx_start=2)

            # 1. test non-streaming
            individuals = parser.parse(stream=False)
            self.assertEqual(
                len(individuals), 1, "Non-stream failed to skip rows"
            )
            self.assertEqual(
                individuals[0].hla_data[0].hla1.allele_string, "A*01:01"
            )

            # 2. test streaming
            individuals_stream = list(parser.parse(stream=True))
            self.assertEqual(
                len(individuals_stream), 1, "Stream failed to skip rows"
            )
            self.assertEqual(
                individuals_stream[0].hla_data[0].hla1.allele_string, "A*01:01"
            )
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_row_idx_start_respected_excel(self):
        """
        Test that row_idx_start is respected in Excel parsing (stream and
        non-stream).
        """
        # create df representing the sheet
        df = pd.DataFrame([
            ["MetaData", "Ignore"],
            ["A_1", "A_2"],
            ["A*01:01", "A*02:01"]
        ])

        with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as f:
            # write without header/index to rely strictly on row positions
            df.to_excel(f.name, index=False, header=False)
            tmp_path = f.name

        try:
            # configure parser to start at Row 2
            parser = HLADataSource(tmp_path, row_idx_start=2)

            # 1. test non-streaming (was failing previously)
            individuals = parser.parse(stream=False)
            self.assertEqual(
                len(individuals), 1, "Non-stream Excel failed to skip rows"
            )
            self.assertEqual(
                individuals[0].hla_data[0].hla1.allele_string, "A*01:01"
            )

            # 2. test streaming
            individuals_stream = list(parser.parse(stream=True))
            self.assertEqual(
                len(individuals_stream), 1, "Stream Excel failed to skip rows"
            )
            self.assertEqual(
                individuals_stream[0].hla_data[0].hla1.allele_string, "A*01:01"
            )
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)


if __name__ == "__main__":
    unittest.main()
