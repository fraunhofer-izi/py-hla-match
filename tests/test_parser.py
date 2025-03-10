import unittest
import os
from py_hla_match.models import Individual
from py_hla_match.parser import HLADataSource


class TestParser(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """Set up file paths for test cases."""
        cls.valid_csv = os.path.join("resources", "hla_test_data.csv")
        cls.valid_excel = os.path.join("resources", "hla_test_data.xlsx")

    def test_parse_valid_csv(self):
        """Test parsing a valid CSV file."""
        parser = HLADataSource(self.valid_csv)
        individuals = parser.parse()
        self.assertIsInstance(individuals, list)
        self.assertEqual(len(individuals), 8)
        self.assertTrue(all(isinstance(ind, Individual) for ind in individuals))
        self.assertEqual(individuals[0].hla_data[4].locus, "DRB1")

    def test_parse_valid_excel(self):
        """Test parsing a valid Excel file."""
        parser = HLADataSource(self.valid_excel)
        individuals = parser.parse()
        self.assertIsInstance(individuals, list)
        self.assertEqual(len(individuals), 8)
        self.assertTrue(all(isinstance(ind, Individual) for ind in individuals))
        self.assertEqual(individuals[0].hla_data[4].locus, "DRB1")


if __name__ == "__main__":
    unittest.main()
