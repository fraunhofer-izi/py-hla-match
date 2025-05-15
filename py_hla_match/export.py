import logging
from itertools import islice

import pandas as pd

from abc import ABC, abstractmethod
from typing import Iterator, Union

from py_hla_match.matching import multi_locus_match
from py_hla_match.models import Individual

logger = logging.getLogger(__name__)


class BaseMatchResult(ABC):

    def __init__(self,
                 source: Union[Iterator[Individual], str],
                 target: Union[Iterator[Individual], str],
                 resolution: str = "basic",
                 chunksize: int = 10000):

        if resolution not in {"basic", "high", "full"}:
            raise ValueError(
                "Resolution must be one of: 'basic', 'high', 'full'"
            )

        self.source = self._load_data(source, chunksize)
        self.target = self._load_data(target, chunksize)
        self.resolution = resolution
        self.chunksize = chunksize
        self.result_file = "match_results.csv"
        self.calculate_result()

    def _load_data(self,
                   data: Union[Iterator[Individual], str],
                   chunksize: int) -> Iterator[pd.DataFrame]:
        """Loads data either from memory or a file."""
        # Data is a file path
        if isinstance(data, str):
            logger.info(f"Loading data from file {data}")
            if data.endswith(".csv"):
                return pd.read_csv(data, chunksize=chunksize)
            elif data.endswith(".xlsx"):
                return pd.read_excel(data, chunksize=chunksize)
            else:
                raise ValueError(
                    "Unsupported file format. Must be either CSV or Excel."
                )
        # Data is a list of Individual objects
        return iter(data)

    @abstractmethod
    def process_chunk(self,
                      source_chunk: Iterator[Individual],
                      target_chunk: Iterator[Individual]) -> pd.DataFrame:
        """Processes a single chunk of data."""
        pass

    def calculate_result(self) -> None:
        """Processes data in chunks and writes to a file."""
        logger.info("Starting result calculation...")
        with open(self.result_file, "w") as f:
            first_chunk = True
            current_chunk = 0
            source_individuals = list(
                islice(
                    self.source,
                    self.chunksize * current_chunk,
                    self.chunksize * (current_chunk + 1)
                )
            )
            target_individuals = list(
                islice(
                    self.target,
                    self.chunksize * current_chunk,
                    self.chunksize * (current_chunk + 1)
                )
            )
            # while there are still chunks to process
            while len(source_individuals) > 0 and len(target_individuals) > 0:
                result = self.process_chunk(
                    source_individuals, target_individuals
                )                
                if first_chunk:
                    result.to_csv(f, index=False)
                    first_chunk = False
                else:
                    result.to_csv(f, index=False, header=False)
                current_chunk += 1
                source_individuals = list(
                    islice(
                        self.source,
                        self.chunksize * current_chunk,
                        self.chunksize * (current_chunk + 1)
                    )
                )
                target_individuals = list(
                    islice(
                        self.target,
                        self.chunksize * current_chunk,
                        self.chunksize * (current_chunk + 1)
                    )
                )
        return pd.read_csv(self.result_file)

    def to_df(self):
        """Loads the final result from the saved CSV."""
        return pd.read_csv(self.result_file)

    def to_csv(self, filename: str) -> None:
        """Saves results to a CSV file."""
        pd.read_csv(self.result_file).to_csv(filename, index=False)

    def to_excel(self, filename: str) -> None:
        """Saves results to an Excel file."""
        pd.read_csv(self.result_file).to_excel(filename, index=False)


class PairwiseMatchResult(BaseMatchResult):
    """
    Match individuals row-wise based on two data sources - > indices of source
    to same indices of target. Use case for this result would be the
    evaluation of an already prepared match resource for several patients, i.e.
    in a study setting.
    """

    def process_chunk(
            self,
            source_chunk: Iterator[Individual],
            target_chunk: Iterator[Individual]
    ) -> pd.DataFrame:
        results = []
        for idx, source_individual in enumerate(source_chunk):
            target_individual = target_chunk[idx]
            match_results = multi_locus_match(
                source_individual, target_individual
            )
            row = {}
            for result in match_results:
                locus = result.patient.locus
                match_level = result.get_match_level_for_resolution(
                    self.resolution
                )
                row[locus] = match_level
            results.append(row)
        df = pd.DataFrame(results)
        return df


class BestMatchResult(BaseMatchResult):
    """
    Find the best match for a single source Individual across a  target
    dataset. Streams through the target in chunks and keeps track of the
    overall best matching donor Individual.
    """

    def __init__(self,
                 source: Individual,
                 target: Union[Iterator[Individual], str],
                 resolution: str = "basic",
                 chunksize: int = 10000):
        # Wrap the single source Individual into a 1‐element iterator
        self._single_source = source
        self.source = iter([source])
        self.target = self._load_data(target, chunksize)
        self.resolution = resolution
        self.chunksize = chunksize
        self.result_file = "best_match_result.csv"
        self.best_match_idx = None
        self.best_match_individual = None
        self.calculate_result()

    def get_best_match(self) -> Individual:
        """
        Returns the best matching donor Individual.
        """
        return self.best_match_individual

    def get_best_match_target_idx(self) -> int:
        """
        Returns the index of the best matching donor Individual in the target
        dataset.
        """
        return self.best_match_idx

    def process_chunk(self) -> pd.DataFrame:
        # not needed, handling logic in calculate result
        return None

    def calculate_result(self) -> None:
        """
        Loads the target data in chunks and processes each chunk to find the
        best match.
        """
        best_score = -1
        source_ind = self._single_source
        with open(self.result_file, "w") as f:  # noqa: F841
            current_chunk = 0
            target_individuals = list(
                islice(
                    self.target,
                    self.chunksize * current_chunk,
                    self.chunksize * (current_chunk + 1)
                )
            )
            # track global index over all chunks
            current_idx = 0
            # loop will stop when new chunk is attempted to be loaded but is
            # empty
            while len(target_individuals) > 0:
                for idx, target_individual in enumerate(target_individuals):
                    match_results = multi_locus_match(
                        source_ind, target_individual
                    )
                    score = sum(
                        result.allele_score for result in match_results
                    )
                    if score > best_score:
                        best_score = score
                        self.best_match_idx = idx
                        self.best_match_individual = target_individual
                current_idx += idx
                current_chunk += 1
                target_individuals = list(
                    islice(
                        self.target,
                        self.chunksize * current_chunk,
                        self.chunksize * (current_chunk + 1)
                    )
                )
