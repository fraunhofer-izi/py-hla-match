import logging
import time
import pandas as pd

from abc import ABC, abstractmethod
from typing import Iterator, Union
from py_hla_match.models import Individual

logger = logging.getLogger(__name__)

class BaseMatchResult(ABC):  
    
    def __init__(self, 
                 source: Union[Iterator[Individual], str], 
                 target: Union[Iterator[Individual], str], 
                 resolution: str = "basic", 
                 chunksize: int = 10000):
        
        if resolution not in {"basic", "high", "full"}:
            raise ValueError("Resolution must be one of: 'basic', 'high', 'full'")
        
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
                raise ValueError("Unsupported file format. Must be either CSV or Excel.")
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
        start_time = time.time()
        with open(self.result_file, "w") as f:
            first_chunk = True
            chunk_count = 0
            for source_chunk in self.source:
                for target_chunk in self.target:
                    chunk_start_time = time.time()
                    result_chunk = self.process_chunk(source_chunk, target_chunk)
                    # flush to file
                    result_chunk.to_csv(f, index=False, header=first_chunk)
                    # Only write header for the first chunk
                    first_chunk = False  
                    chunk_count += 1
                    logger.info(f"Processed chunk {chunk_count}: {len(source_chunk)} source, {len(target_chunk)} target records in {time.time() - chunk_start_time:.2f} seconds")
        logger.info(f"Completed result calculation in {time.time() - start_time:.2f} seconds. Results saved to {self.result_file}")

    def to_df(self) -> pd.DataFrame:
        """Loads the final result from the saved CSV."""
        return pd.read_csv(self.result_file)

    def to_csv(self, filename: str) -> None:
        """Saves results to a CSV file."""
        pd.read_csv(self.result_file).to_csv(filename, index=False)

    def to_excel(self, filename: str) -> None:
        """Saves results to an Excel file."""
        pd.read_csv(self.result_file).to_excel(filename, index=False)

        
class PairwiseMatchResult(BaseMatchResult):
    
    def process_chunk(self, source_chunk: pd.DataFrame, target_chunk: pd.DataFrame) -> pd.DataFrame:
        results = []
        for _, source_row in source_chunk.iterrows():
            for _, target_row in target_chunk.iterrows():
                score = self._compare(source_row, target_row)  # Implement your matching logic
                results.append({"source_id": source_row["id"], "target_id": target_row["id"], "score": score})
        return pd.DataFrame(results)

    def _compare(self, source: pd.Series, target: pd.Series) -> float:
        """Placeholder comparison logic (to be implemented)."""
        return 1.0  # Dummy score


class BestMatchResult(BaseMatchResult):
    
    def __init__(self, source: Union[Iterator[Individual], str], target: Union[Iterator[Individual], str], 
                 resolution: str = "basic", n_candidates: int = 1, chunksize: int = 10000):
        super().__init__(source, target, resolution, chunksize)
        self.n_candidates = n_candidates

    def process_chunk(self, source_chunk: pd.DataFrame, target_chunk: pd.DataFrame) -> pd.DataFrame:
        """Finds the best `n_candidates` matches per source individual."""
        results = []
        for _, source_row in source_chunk.iterrows():
            scores = []
            for _, target_row in target_chunk.iterrows():
                score = self._compare(source_row, target_row)
                scores.append({"source_id": source_row["id"], "target_id": target_row["id"], "score": score})
            
            top_matches = sorted(scores, key=lambda x: x["score"], reverse=True)[:self.n_candidates]
            results.extend(top_matches)
        
        return pd.DataFrame(results)

    def _compare(self, source: pd.Series, target: pd.Series) -> float:
        """Placeholder for scoring logic."""
        return 1.0  # Dummy score

