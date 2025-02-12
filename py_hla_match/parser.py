import logging

from py_hla_match.exceptions import MalformedHLAStringError, MalformedHLADataSourceError
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
        individuals: [Individual] = []
        # slice the dataframe if start and end indices were given
        if self.col_idx_start and self.col_idx_stop:
            df = df.iloc[:, self.col_idx_start:self.col_idx_stop]
        for idx, row in df.iterrows():
            hla_pairs: [HLAPair] = []
            individual_hla_objects: HLA = []
            # first: parse all available HLA strings in the row into HLA objects
            for hla_string in row:
                try:
                    hla = HLA(hla_string)
                    individual_hla_objects.append(hla)
                except MalformedHLAStringError as e:
                    logger.error(f'Encountered malformed HLA String {hla_string} in row {idx}. Skipping Allele.')
                    continue
            # now: Match HLA pairs based on locus
            for hla_object in individual_hla_objects:
                # check if we already have a pair for the locus
                for pair in hla_pairs:
                    if pair.hla1.locus == hla_object.locus:
                        # edge case: we already got two alleles for the locus, but there is a third one in the data
                        if pair.hla2 is not None:
                            raise MalformedHLADataSourceError(
                                f"Encountered third allele for locus {hla_object.locus} in row {idx}.")
                        pair.hla2 = hla_object
                        break
                else:
                    # if we don't have a pair for the locus, create a new pair
                    hla_pairs.append(HLAPair(hla1=hla_object, hla2=None))
                # matched -> remove from individual_hla_objects
                individual_hla_objects.remove(hla_object)
            # edge case: check for unpaired allele
            for hla_pair in hla_pairs:
                if hla_pair.hla2 is None:
                    logger.warning(f"Unpaired allele {hla_pair.hla1.allele_string} in row {idx}.")
            individuals.append(Individual(hla_data=hla_pairs))
            logger.info(f"Successfully parsed row {idx}. Added {len(hla_pairs)} HLA pairs to individual.")
        return individuals
