import logging
from itertools import zip_longest

import pandas as pd

from typing import List
from pandas import DataFrame

from py_hla_match.matching import MatchResult, multi_locus_match
from py_hla_match.parser import HLADataSource

logger = logging.getLogger(__name__)


class PairwiseMatchResult:
    """
    Match individuals row-wise based on two data sources - > indices of source
    to same indices of target. Will store the results in a csv file.

    :param source: HLADataSource for the source dataset
    :param target: HLADataSource for the target dataset
    :param storage_filename: Name of the file to store the results
    :param resolution: Resolution of the match results, can be 'basic', 'high', or 'full'
    :param stream: If True, results will be streamed and not stored in memory
    :param chunk_size: Size of the chunks to read from the file (if streaming)

    :raises ValueError: If resolution is not one of 'basic', 'high', or 'full'
    """

    def __init__(self,
                 source: HLADataSource,
                 target: HLADataSource,
                 storage_filename: str = "match_results.csv",
                 resolution: str = "basic",
                 stream: bool = False,
                 chunk_size: int = 10000):

        if resolution not in {"basic", "high", "full"}:
            raise ValueError(
                "Resolution must be one of: 'basic', 'high', 'full'"
            )

        self.source = source
        self.target = target
        self.resolution = resolution
        self.stream = stream
        self.chunk_size = chunk_size
        self.result_file = storage_filename
        self.result = None  # Placeholder for the result DataFrame
        self.calculate_result()


    def to_df(self) -> DataFrame:
        """
        Converts the match results to a pandas DataFrame.

        :return: DataFrame containing the match results
        """
        if self.stream:
            raise ValueError("Cannot convert to DataFrame when streaming is enabled.")
        return self.result

    def calculate_result(self) -> None:
        """
        Matches individuals from source and target datasets row-wise.
        Assumes that both datasets are aligned by index.
        Processes data in chunks and periodically flushes results to the output file.
        """
        logger.info("Starting pairwise match result calculation...")

        # Parse source and target data
        source_data = self.source.parse(stream=self.stream, chunk_size=self.chunk_size)
        target_data = self.target.parse(stream=self.stream, chunk_size=self.chunk_size)

        # Handle non-streaming case by converting lists to single-chunk iterators
        if not self.stream:
            # Flatten the lists into iterables of Individual objects
            source_data = iter(source_data)
            target_data = iter(target_data)

        all_loci = set()
        # Keep track of headers in memory
        current_headers = None
        is_first_chunk = True

        for idx, (source_ind, target_ind) in enumerate(zip_longest(source_data, target_data, fillvalue=None)):

            if source_ind is None:
                raise ValueError(f"source_data exhausted before target_data at index {idx}")
            elif target_ind is None:
                raise ValueError(f"target_data exhausted before source_data at index {idx}")

            results = []

            # Process individual objects directly
            match_results: List[MatchResult] = multi_locus_match(source_ind, target_ind)
            row = {}
            for result in match_results:
                locus = result.patient.locus
                row[locus] = result.get_match_level_for_resolution(self.resolution)
                all_loci.add(locus)
            results.append(row)

            # Create a DataFrame for the current chunk
            chunk_df = pd.DataFrame(results)
            chunk_df = chunk_df.reindex(columns=sorted(all_loci))  # Ensure all loci are included

            # Write to the results file
            if self.stream:
                new_headers = sorted(all_loci)
                if is_first_chunk:
                    # Write header for the first chunk
                    chunk_df.to_csv(self.result_file, index=False, mode='w')
                    current_headers = new_headers
                    is_first_chunk = False
                else:
                    # Check if the header needs to be updated
                    if current_headers != new_headers:
                        # Rewrite the file with updated headers
                        existing_data = pd.read_csv(self.result_file)
                        existing_data = existing_data.reindex(columns=new_headers, fill_value=None)
                        existing_data.to_csv(self.result_file, index=False, mode='w')
                        chunk_df.to_csv(self.result_file, index=False, mode='a', header=False)
                        current_headers = new_headers
                    else:
                        # Append without writing the header for subsequent chunks
                        chunk_df.to_csv(self.result_file, index=False, mode='a', header=False)
            else:
                # If streaming is disabled, accumulate results in memory
                if self.result is None:
                    self.result = chunk_df
                else:
                    self.result = pd.concat([self.result, chunk_df], ignore_index=True)

        # Write the accumulated results to the file in non-streaming mode
        if not self.stream and self.result is not None:
            self.result.to_csv(self.result_file, index=False)

        logger.info("Pairwise match result calculation completed.")
