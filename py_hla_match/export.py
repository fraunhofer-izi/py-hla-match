import pandas as pd

from abc import ABC, abstractmethod
from typing import List

from py_hla_match.models import Individual

class BaseMatchResult(ABC):  
    def __init__(self, source: List[Individual], target: List[Individual], resolution: str = "basic"):
        if resolution not in {"basic", "high", "full"}:
            raise ValueError("Resolution must be one of: 'basic', 'high', 'full'")
        self.source = source
        self.target = target
        self.resolution = resolution
        self.result = self.calculate_result()
        
    @abstractmethod
    def calculate_result(self) -> pd.DataFrame:
        """Calculates the match result. Must be implemented by subclasses."""
        pass

    def to_df(self) -> pd.DataFrame:
        """Exports the result to a Pandas DataFrame."""
        return pd.DataFrame(self.result)

    def to_csv(self, filename: str) -> None:
        """Exports the result to a CSV file."""
        self.to_df().to_csv(filename, index=False)

    def to_excel(self, filename: str) -> None:
        """Exports the result to an Excel file."""
        self.to_df().to_excel(filename, index=False)
        
class PairwiseMatchResult(BaseMatchResult):
    def __init__(self, source: List[Individual], target: List[Individual], resolution: str = "basic"):
        super().__init__(source, target, resolution)

class BestMatchResult(BaseMatchResult):
    def __init__(self, source: List[Individual], target: List[Individual], resolution: str = "basic"):
        super().__init__(source, target, resolution)
