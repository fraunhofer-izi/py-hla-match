Quickstart
==========

This quickstart uses the artificial CSVs under the demo folder and avoids any
real or sensitive data.

Run a basic pairwise match
--------------------------

Use the synthetic patient and donor CSVs and write results to a new file:

.. code-block:: python

   from py_hla_match.parser import HLADataSource
   from py_hla_match.export import PairwiseMatch

   data_path = "py_hla_match/demo/data/random_data/synthetic_patients.csv"
   donor_path = "py_hla_match/demo/data/random_data/synthetic_donors.csv"
   output_path = "py_hla_match/demo/data/random_data/match_results.csv"

   src = HLADataSource(
       data_path,
       col_idx_start=1,
       col_idx_stop=13,
       row_idx_start=1,
   )

   tgt = HLADataSource(
       donor_path,
       col_idx_start=1,
       col_idx_stop=13,
       row_idx_start=1,
   )

   matcher = PairwiseMatch(
       source=src,
       target=tgt,
       storage_filename=output_path,
       resolution="high",
   )

   matcher.run()

Inspect raw allele-level results
--------------------------------

Convert raw match levels to a DataFrame and write to CSV:

.. code-block:: python

   raw_output_path = "py_hla_match/demo/data/random_data/match_results_raw.csv"
   matcher.raw_to_df().to_csv(raw_output_path, index=False)
