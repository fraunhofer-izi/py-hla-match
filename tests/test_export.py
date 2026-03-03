import os
import unittest
import tempfile
import shutil
import pandas as pd
import logging  # noqa: F401

from py_hla_match.export import PairwiseMatch
from py_hla_match.parser import HLADataSource


def _generate_tmp_file_source_and_target(
        source_df: pd.DataFrame,
        target_df: pd.DataFrame
) -> tuple[HLADataSource, HLADataSource]:
    """
    Generate temporary files for source and target HLA data.
    """
    with tempfile.NamedTemporaryFile(
        delete=False, suffix=".csv", mode='w'
    ) as sf:
        source_df.to_csv(sf.name, index=False)
        source_path = sf.name
    with tempfile.NamedTemporaryFile(
        delete=False, suffix=".csv", mode='w'
    ) as tf:
        target_df.to_csv(tf.name, index=False)
        target_path = tf.name

    # NOTE.: HLADataSource needs path
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
        if os.path.exists(cls.temp_dir):
            shutil.rmtree(cls.temp_dir)

    def tearDown(self):
        """Clean up generated files."""
        # cleanup result files created during tests
        for filename in os.listdir(self.temp_dir):
            file_path = os.path.join(self.temp_dir, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)

    def test_valid_csv_no_streaming_defaults(self):
        """
        Test parsing valid CSV without streaming using default settings.
        Should include all canonical loci.
        """
        source = HLADataSource(self.valid_csv)
        target = HLADataSource(self.valid_csv)
        output_file = os.path.join(self.temp_dir, "match_defaults.csv")

        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=output_file,
            stream=False
        )
        pairwise_match.run()

        self.assertTrue(os.path.exists(output_file))

        # check content
        result_df = pairwise_match.to_df()
        self.assertIsInstance(result_df, pd.DataFrame)

        # validation: should have at least A, B, C, DRB1...
        # just checking existence of columns for basic output
        self.assertIn("A_1", result_df.columns)
        self.assertIn("A_2", result_df.columns)
        self.assertIn("DRB1_1", result_df.columns)

        # should not have detail columns by default
        self.assertNotIn("A_ard_1", result_df.columns)
        self.assertNotIn("A_mol_1", result_df.columns)

    def test_feature_flags_include_details(self):
        """
        Test that enabling feature flags adds the correct columns.
        """
        source = HLADataSource(self.valid_csv)
        target = HLADataSource(self.valid_csv)
        output_file = os.path.join(self.temp_dir, "match_details.csv")

        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=output_file,
            stream=False,
            include_ard_details=True,
            include_molecular_details=True,
            include_homozygosity=True
        )
        pairwise_match.run()

        result_df = pairwise_match.to_df()

        # check ard columns
        self.assertIn("A_ard_1", result_df.columns)
        self.assertIn("A_ard_cert_1", result_df.columns)

        # check molecular columns
        self.assertIn("A_mol_1", result_df.columns)

        # check homozygosity
        self.assertIn("A_homozygous_patient", result_df.columns)

    def test_explicit_loci_filtering(self):
        """
        Test that providing 'loci' list filters the output.
        """
        source = HLADataSource(self.valid_csv)
        target = HLADataSource(self.valid_csv)
        output_file = os.path.join(self.temp_dir, "match_filtered.csv")

        # request 'A' only, even though file has A, B, C, DRB1
        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=output_file,
            loci=["A"],
            stream=False
        )
        pairwise_match.run()

        result_df = pairwise_match.to_df()

        # must only have A
        self.assertIn("A_1", result_df.columns)

        # and not anything else
        self.assertNotIn("B_1", result_df.columns)
        self.assertNotIn("C_1", result_df.columns)
        self.assertNotIn("DRB1_1", result_df.columns)

    def test_overwrite_protection(self):
        """
        Test that overwrite=False raises FileExistsError.
        """
        output_file = os.path.join(self.temp_dir, "existing.csv")
        # create dummy file
        with open(output_file, 'w') as f:
            f.write("dummy")

        source = HLADataSource(self.valid_csv)
        target = HLADataSource(self.valid_csv)

        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=output_file,
            overwrite=False
        )

        with self.assertRaises(FileExistsError):
            pairwise_match.run()

        # should work if overwrite=True
        pairwise_match_overwrite = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=output_file,
            overwrite=True
        )
        try:
            pairwise_match_overwrite.run()
        except FileExistsError:
            self.fail(
                "run() raised FileExistsError unexpectedly with overwrite=True"
            )

    def test_valid_csv_streaming(self):
        """Test parsing a valid CSV file with streaming."""
        source = HLADataSource(self.valid_csv)
        target = HLADataSource(self.valid_csv)
        output_file = os.path.join(self.temp_dir, "match_streaming.csv")

        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=output_file,
            stream=True,
            chunk_size=2,
            overwrite=True
        )
        pairwise_match.run()

        # assert file exists
        self.assertTrue(os.path.exists(output_file))

        # assert exception on to_df
        with self.assertRaises(RuntimeError):
            pairwise_match.to_df()

        # verify content via manual read
        df = pd.read_csv(output_file)
        self.assertGreater(len(df), 0)
        self.assertIn("A_1", df.columns)

    def test_length_mismatch_raises(self):
        """
        PairwiseMatchResult must raise a ValueError on length mismatch.
        """
        source_df = pd.DataFrame(
            {"A1": ["A*01:01", "A*01:01"],
             "A2": ["A*02:01", "A*02:01"]}
        )  # 2 rows
        target_df = pd.DataFrame(
            {"A1": ["A*01:01", "A*01:01", "A*01:01"],
             "A2": ["A*02:01", "A*02:01", "A*02:01"]}
        )  # 3 rows

        source, target = _generate_tmp_file_source_and_target(
            source_df, target_df
        )

        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(
                self.temp_dir, "len_mismatch.csv"
            ),
            stream=False
        )

        with self.assertRaises(ValueError):
            pairwise_match.run()

    def test_unexpected_locus_handling(self):
        """
        Tests behavior when a locus appears in data that isn't expected.
        """
        # has A and C (unexpected if we only look for A and B)
        source_df = pd.DataFrame({
             "A1": ["A*01:01"], "A2": ["A*02:01"],
             "C1": ["C*01:02"], "C2": ["C*07:02"]
        })
        # Target has A and B
        target_df = pd.DataFrame({
             "A1": ["A*01:01"], "A2": ["A*02:01"],
             "B1": ["B*07:02"], "B2": ["B*44:02"]
        })

        source, target = _generate_tmp_file_source_and_target(
            source_df, target_df
        )

        output_file = os.path.join(self.temp_dir, "mixed_loci.csv")

        # 1: default (all of loci)
        # should capture A, B, and C and all other loci
        pm_default = PairwiseMatch(source, target, output_file, overwrite=True)
        pm_default.run()
        df_default = pm_default.to_df()

        self.assertIn("A_1", df_default.columns)
        self.assertIn("B_1", df_default.columns)
        self.assertIn("C_1", df_default.columns)

        # 2: filter (A and B only)
        # should capture A and B, drop C
        pm_filter = PairwiseMatch(
            source, target, output_file,
            loci=["A", "B"], overwrite=True
        )
        pm_filter.run()
        df_filter = pm_filter.to_df()

        self.assertIn("A_1", df_filter.columns)
        self.assertIn("B_1", df_filter.columns)
        self.assertNotIn("C_1", df_filter.columns)

    def test_invalid_csv_no_streaming(self):
        """Test parsing an invalid CSV file without streaming."""
        source = HLADataSource(self.invalid_csv)
        target = HLADataSource(self.invalid_csv)

        pairwise_match = PairwiseMatch(
            source=source,
            target=target,
            storage_filename=os.path.join(
                self.temp_dir, "match_invalid.csv"
            ),
            stream=False
        )

        # check that it logs errors but proceeds
        with self.assertLogs("py_hla_match.parser", level="ERROR") as cm:
            pairwise_match.run()
            self.assertTrue(
                any("malformed HLA String" in r.message for r in cm.records)
            )

        result_df = pairwise_match.to_df()
        self.assertGreater(len(result_df), 0)
