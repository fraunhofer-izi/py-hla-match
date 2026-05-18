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

   data_path = "demo/data/random_data/synthetic_patients.csv"
   donor_path = "demo/data/random_data/synthetic_donors.csv"
   output_path = "demo/data/random_data/match_results.csv"

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
       include_ard_details=True,
       include_molecular_details=True,
       include_dpb1_tce=False,
       include_homozygosity=False,
       overwrite=True,
   )

   matcher.run()

Explore results
---------------

Results are written to ``output_path`` during ``matcher.run()``. To inspect
in-memory:

.. code-block:: python

   df = matcher.to_df()
   print(df.head())