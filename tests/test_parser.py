import logging
import unittest
import os
from py_hla_match.models import Individual
from py_hla_match.parser import HLADataSource


class TestParser(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up file paths for test cases."""
        TEST_DIR_PATH = os.path.dirname(os.path.realpath(__file__))
        cls.valid_csv = os.path.join(TEST_DIR_PATH, "resources", "hla_test_data.csv")
        cls.valid_excel = os.path.join(TEST_DIR_PATH, "resources", "hla_test_data.xlsx")
        cls.invalid_csv = os.path.join(TEST_DIR_PATH, "resources", "hla_test_data_malformed.csv")

    def test_parse_valid_csv(self):
        """Test parsing a valid CSV file."""
        parser = HLADataSource(self.valid_csv)
        individuals = parser.parse()
        self.assertIsInstance(individuals, list)
        self.assertEqual(len(individuals), 8)
        self.assertTrue(all(isinstance(ind, Individual) for ind in individuals))
        self.assertEqual("DRB1", individuals[0].hla_data[4].locus)

    def test_parse_malformed_hla_strings_get_logged(self):
        """Test parsing an invalid CSV file logs MalformedHLAStringError twice."""
        parser = HLADataSource(self.invalid_csv)
        with self.assertLogs("py_hla_match.parser", level="ERROR") as log_context:
            parser.parse()
        error_logs = [record for record in log_context.output if "Encountered malformed HLA String" in record]
        self.assertEqual(len(error_logs), 3, "Expected exactly two malformed HLA string log entries")

    def test_parse_valid_excel(self):
        """Test parsing a valid Excel file."""
        parser = HLADataSource(self.valid_excel)
        individuals = parser.parse()
        self.assertIsInstance(individuals, list)
        self.assertEqual(len(individuals), 8)
        self.assertTrue(all(isinstance(ind, Individual) for ind in individuals))
        self.assertEqual("DRB1", individuals[0].hla_data[4].locus)


if __name__ == "__main__":
    unittest.main()
