import logging
import os
from itertools import zip_longest

import pandas as pd
import csv

from typing import List, Dict, Any, Optional, Iterable

from py_hla_match.config import get_config, LOCUS_ALIAS_MAP
from py_hla_match.matching import MatchResult, multi_locus_match
from py_hla_match.parser import HLADataSource

logger = logging.getLogger(__name__)


def scan_loci(source: HLADataSource, chunk_size: int = 10000) -> List[str]:
    """
    Utility function to scan hla data source and identify all loci present.

    :param source: HLADataSource to scan
    :param chunk_size: Size of the chunks to read from the file
    :return: Sorted list of unique loci detected in the data source
    """
    detected = set()
    # we stream through the file to avoid memory overhead during scan
    for ind in source.parse(stream=True, chunk_size=chunk_size):
        for pair in ind.hla_data:
            detected.add(pair.locus)
    return sorted(list(detected))


def _flatten_match_results(
    results: List[MatchResult],
    include_ard: bool = False,
    include_molecular: bool = False,
    include_homozygosity: bool = False,
    include_dpb1_tce: bool = False,
) -> Dict[str, Any]:
    """
    Transforms a list of MatchResult objects into a flat dictionary.

    Output keys use raw Enum names (e.g., 'A_1': 'ARD_MATCH').

    :param results: List of MatchResult objects for one pair
    :param include_ard: If True, adds ARD match levels and certainties^
    :param include_molecular: If True, adds molecular levels and certainties
    :param include_homozygosity: If True, adds patient homozygosity status
    :param include_dpb1_tce: If True, adds DPB1 TCE status (if computed)
    :return: Dictionary representing a single flattened row
    """
    row = {}
    for res in results:
        locus = res.patient.locus

        # identify target keys (loci)
        target_loci = [locus]
        if locus in LOCUS_ALIAS_MAP:
            target_loci.extend(LOCUS_ALIAS_MAP[locus])

        # populate dictionary for all targets
        for target_locus in target_loci:
            # base Levels
            row[f"{target_locus}_1"] = res.allele_match_levels[0].name
            row[f"{target_locus}_2"] = res.allele_match_levels[1].name

            # pptional details
            if include_ard:
                row[f"{target_locus}_ard_1"] = \
                    res.ard_match_levels[0].name
                row[f"{target_locus}_ard_cert_1"] = \
                    res.ard_match_certainties[0].name
                row[f"{target_locus}_ard_2"] = \
                    res.ard_match_levels[1].name
                row[f"{target_locus}_ard_cert_2"] = \
                    res.ard_match_certainties[1].name

            if include_molecular:
                row[f"{target_locus}_mol_1"] = \
                    res.molecular_match_levels[0].name
                row[f"{target_locus}_mol_cert_1"] = \
                    res.molecular_match_certainties[0].name
                row[f"{target_locus}_mol_2"] = \
                    res.molecular_match_levels[1].name
                row[f"{target_locus}_mol_cert_2"] = \
                    res.molecular_match_certainties[1].name

            if include_homozygosity and res.is_homozygous_patient is not None:
                row[f"{target_locus}_homozygous_patient"] = \
                    res.is_homozygous_patient

        # DPB1 TCE (Only for DPB1)
        if include_dpb1_tce and locus == "DPB1" and res.dpb1_tce_status:
            row["DPB1_tce_status"] = res.dpb1_tce_status.name

    return row


class PairwiseMatch:
    """
    Match individuals row-wise based on two data sources - > indices of source
    to same indices of target. Will store the results in a csv file.

    :param source: HLADataSource for the source dataset
    :param target: HLADataSource for the target dataset
    :param storage_filename: Name of the file to store the results
    :param loci: Optional iterable of specific loci to export. If None,
        defaults to all supported loci
    :param include_ard_details: If True, include ARD refinement columns
    :param include_molecular_details: If True, include molecular refinement
        columns
    :param include_homozygosity: If True, include homozygosity boolean
    :param include_dpb1_tce: If True, include DPB1 TCE status column
    :param stream: If True, results will be streamed and not stored in
        memory
    :param chunk_size: Size of the chunks to read from the file
    :param overwrite: If True, allow overwriting existing output files

    :raises ValueError: If resolution is not one of 'basic', 'high', or 'full'
    """
    def __init__(
        self,
        source: HLADataSource,
        target: HLADataSource,
        storage_filename: str = "match_results.csv",
        # Config
        loci: Optional[Iterable[str]] = None,
        # Feature Flags
        include_ard_details: bool = False,
        include_molecular_details: bool = False,
        include_homozygosity: bool = False,
        include_dpb1_tce: bool = False,
        # Controls
        stream: bool = False,
        chunk_size: int = 10000,
        overwrite: bool = False
    ):
        self.source = source
        self.target = target
        self.storage_filename = storage_filename

        # define loci upfront
        if loci:
            # respect user order exactly
            self.loci = list(loci)
        else:
            self.loci = sorted(list(get_config().effective_valid_loci))

        self.include_ard = include_ard_details
        self.include_mol = include_molecular_details
        self.include_hom = include_homozygosity
        self.include_tce = include_dpb1_tce

        self.stream = stream
        self.chunk_size = chunk_size
        self.overwrite = overwrite

        # State (Only used if stream=False)
        self.raw_results: List[List[MatchResult]] = []
        self._df_buffer: List[dict] = []

        self._headers = self._build_headers()

    def _build_headers(self) -> List[str]:
        """
        Builds the CSV header row based on configured loci and flags.

        :return: List of column names
        """
        headers = ['pair_index']
        for locus in self.loci:
            headers.extend([f"{locus}_1", f"{locus}_2"])
            if self.include_ard:
                headers.extend([f"{locus}_ard_1", f"{locus}_ard_cert_1",
                                f"{locus}_ard_2", f"{locus}_ard_cert_2"])
            if self.include_mol:
                headers.extend([f"{locus}_mol_1", f"{locus}_mol_cert_1",
                                f"{locus}_mol_2", f"{locus}_mol_cert_2"])
            if self.include_hom:
                headers.append(f"{locus}_homozygous_patient")

        if self.include_tce and "DPB1" in self.loci:
            headers.append("DPB1_tce_status")
        return headers

    def run(self) -> None:
        """
        Executes matching pipeline.

        Matches individuals from source and target datasets row-wise.
        Processes data in chunks (if streamed) or in memory.

        :raises FileExistsError: If output file exists and overwrite is False
        :raises ValueError: If input datasets have mismatched lengths
        """
        logger.info("Starting match calculation...")

        # check directory and file existence
        result_dir = os.path.dirname(self.storage_filename)
        if result_dir and not os.path.exists(result_dir):
            raise FileNotFoundError(
                f"Directory '{result_dir}' does not exist."
            )

        if not self.overwrite and os.path.exists(self.storage_filename):
            raise FileExistsError(f"File '{self.storage_filename}' exists.")

        # open file
        with open(self.storage_filename, 'w', newline='') as f:
            # extrasaction='ignore' drops loci not in self.loci
            writer = csv.DictWriter(
                f, fieldnames=self._headers, extrasaction='ignore'
            )
            writer.writeheader()

            # seset iterators
            src_iter = self.source.parse(
                stream=self.stream, chunk_size=self.chunk_size
            )
            tgt_iter = self.target.parse(
                stream=self.stream, chunk_size=self.chunk_size
            )

            if not self.stream:
                src_iter, tgt_iter = iter(src_iter), iter(tgt_iter)

            # process
            for idx, (src, tgt) in enumerate(zip_longest(src_iter, tgt_iter)):
                if src is None or tgt is None:
                    raise ValueError(f"Length mismatch at index {idx}")

                results = multi_locus_match(src, tgt)

                # dpb1 tce if requested
                if self.include_tce:
                    for res in results:
                        if res.patient.locus == "DPB1":
                            res.get_dpb1_tce_status()

                # flatten
                row = _flatten_match_results(
                    results,
                    self.include_ard, self.include_mol,
                    self.include_hom, self.include_tce
                )
                row['pair_index'] = idx

                # write
                writer.writerow(row)

                # memory storage
                if not self.stream:
                    self.raw_results.append(results)
                    self._df_buffer.append(row)

        logger.info("Matching completed.")

    def to_df(self) -> pd.DataFrame:
        """
        Returns a DataFrame of the results.

        Only available if stream=False.

        :return: pandas DataFrame containing the match results
        :raises RuntimeError: If streaming is enabled
        """
        if self.stream:
            raise RuntimeError("to_df() not available in streaming mode.")
        if not self._df_buffer:
            return pd.DataFrame(columns=self._headers)
        return pd.DataFrame(self._df_buffer).reindex(columns=self._headers)
