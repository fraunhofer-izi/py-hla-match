import logging
from collections import defaultdict

from py_hla_match.exceptions import MalformedHLAStringError, \
    MalformedHLADataSourceError
from py_hla_match.hla import HLA
from py_hla_match.models import Individual, HLAPair

import pandas as pd

logger = logging.getLogger(__name__)


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

    def parse(self) -> list[Individual]:
        """
        Parse HLA data from an excel or csv file.
        """
        if self.source_path.endswith('.xlsx'):
            return self._parse_excel()
        if self.source_path.endswith('.csv'):
            return self._parse_csv()
        raise ValueError("Unsupported file format.")

    def _parse_excel(self) -> list[Individual]:
        """
        Parse HLA data from an excel file.
        """
        df = pd.read_excel(self.source_path)
        return self._parse_dataframe(df)

    def _parse_csv(self) -> list[Individual]:
        """
        Parse HLA data from a csv file.
        """
        df = pd.read_csv(self.source_path)
        return self._parse_dataframe(df)

    def _parse_dataframe(self, df: pd.DataFrame) -> list[Individual]:
        """
        Parse HLA data from a pandas DataFrame.
        """
        individuals: list[Individual] = []
        # slice the dataframe if start and end indices were given
        if self.col_idx_start and self.col_idx_stop:
            df = df.iloc[:, self.col_idx_start:self.col_idx_stop]
        for idx, row in df.iterrows():
            hla_pairs: list[HLAPair] = []
            # Map of locus to HLA objects
            locus_map = defaultdict(list)
            individual_hla_objects: list[HLA] = []
            # first: parse all available HLA strings in the row into HLA
            # objects
            for hla_string in row:
                try:
                    hla = HLA(hla_string)
                    individual_hla_objects.append(hla)
                    locus_map[hla.locus].append(hla)
                except MalformedHLAStringError as e:  # NOQA
                    logger.error(
                        f'Encountered malformed HLA String {hla_string} in row'
                        f' {idx}. Skipping Allele.'
                    )
                    continue
            # now: Match HLA pairs based on locus
            for locus, alleles in locus_map.items():
                # edge case: more than two alleles found for a locus
                if len(alleles) > 2:
                    raise MalformedHLADataSourceError(
                        f"Encountered third allele for locus {locus} in row"
                        f"{idx}."
                    )
                # only create a pair if exactly two alleles exist
                if len(alleles) == 2:
                    hla_pairs.append(HLAPair(hla1=alleles[0], hla2=alleles[1]))
                else:
                    logger.warning(
                        f"Unpaired allele {alleles[0].allele_string} in row"
                        f" {idx}."
                    )
            individuals.append(Individual(hla_data=hla_pairs))
            logger.info(
                f"Successfully parsed row {idx}. Added {len(hla_pairs)} HLA"
                f" pairs to individual."
            )
        return individuals
