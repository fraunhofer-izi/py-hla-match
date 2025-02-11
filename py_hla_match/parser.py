from py_hla_match.hla import HLA
from py_hla_match.models import Individual, HLAPair

import pandas as pd


class HLADataSource:
    """
    Data source for HLA data. Parses HLA data from an excel or csv file.
    """

    def __init__(self, source_path: str,
                 col_idx_start: int = None,
                 col_idx_stop: int = None) -> None:
        """
        Initialize the HLADataSource.

        :param source_path: Path to the excel or csv file
        :param col_idx_start: Column index to start parsing from
        :param col_idx_stop: Column index to stop parsing at
        """
        self.source_path = source_path
        self.col_idx_start = col_idx_start
        self.col_idx_stop = col_idx_stop

    def parse(self) -> [Individual]:
        """
        Parse HLA data from an excel or csv file.
        """
        if self.source_path.endswith('.xlsx'):
            return self._parse_excel()
        if self.source_path.endswith('.csv'):
            return self._parse_csv()
        raise ValueError("Unsupported file format.")

    def _parse_excel(self) -> [Individual]:
        """
        Parse HLA data from an excel file.
        """
        df = pd.read_excel(self.source_path)
        return self._parse_dataframe(df)

    def _parse_csv(self, source_path: str) -> [Individual]:
        """
        Parse HLA data from a csv file.
        """
        df = pd.read_csv(self.source_path)
        return self._parse_dataframe(df)

    def _parse_dataframe(self, df: pd.DataFrame) -> [Individual]:
        """
        Parse HLA data from a pandas DataFrame.
        """
        individuals = []
        # slice the dataframe if start and end indices were given
        if self.col_idx_start and self.col_idx_stop:
            df = df.iloc[:, self.col_idx_start:self.col_idx_stop]
        for idx, row in df.iterrows():
            hla_data = []
            for col in df.columns:
                hla_data.append(HLAPair(HLA(row[col])))
            individuals.append(Individual(hla_data))
        return individuals
