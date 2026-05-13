import csv
import io
import os
import random
from typing import Dict, List


class SyntheticDataGenerator:
    """
    Generate synthetic HLA typing datasets from haplotype frequencies.

    An individual's full HLA type is assembled by drawing two haplotypes
    from the provided frequency distribution (one per chromosome),
    yielding realistic heterozygous and occasionally homozygous typings
    across six loci.

    :param haplotype_freq_data: Raw text content of a haplotype
        frequency file. Expected CSV format with columns for haplotype
        string (index 1) and frequency float (index 3).
    :type haplotype_freq_data: str

    :raises ValueError: If no valid frequency data can be parsed from
        *haplotype_freq_data*.
    """

    def __init__(self, haplotype_freq_data: str):
        self.haplotypes: List[str] = []
        self.frequencies: List[float] = []
        self._parse_haplotype_data(haplotype_freq_data)

        # define loci in the order they appear in haplotype string
        self.loci = ["A", "B", "C", "DRB1", "DQB1", "DPB1"]

    def _parse_haplotype_data(self, data: str):
        """
        Parse raw haplotype frequency text into internal lists.

        Reads a CSV-formatted string, extracts haplotype identifiers and
        their associated frequencies, filters entries outside the
        ``[0.0, 1.0]`` range, and normalises the remaining frequencies
        to sum to 1.0.

        :param data: Raw CSV text with header row. Column at index 1
            contains the hyphen-delimited haplotype string, column at
            index 3 contains the frequency as a float.
        :type data: str

        :raises ValueError: If no rows with valid frequencies remain
            after filtering.
        """
        # use io.StringIO to treat the string data as a file
        data_file = io.StringIO(data)
        reader = csv.reader(data_file)

        # skip header
        next(reader, None)

        haplotypes_raw = []
        frequencies_raw = []

        for row in reader:
            if not row:
                continue
            try:
                haplotype_str = row[1]
                # some frequencies in the sample are > 1.0, which is not
                # possible. We will collect all valid floats and then
                # normalize them.
                frequency = float(row[3])
                if 0.0 <= frequency <= 1.0:  # Simple validation
                    haplotypes_raw.append(haplotype_str)
                    frequencies_raw.append(frequency)
            except (IndexError, ValueError):
                # skip malformed rows
                pass

        # normalize frequencies so they sum to 1, making them a valid
        # probability distribution
        total_freq = sum(frequencies_raw)
        if total_freq > 0:
            self.haplotypes = haplotypes_raw
            self.frequencies = [f / total_freq for f in frequencies_raw]
        else:
            raise ValueError("No valid frequency data could be parsed.")

    def _generate_individual(self) -> Dict[str, str]:
        """
        Generate a full HLA type for a single individual.

        Draws two haplotypes from the frequency distribution (with
        replacement) and assigns alleles to each locus for both
        chromosomes.

        :returns: HLA type keyed by ``"{locus}*1"`` and ``"{locus}*2"``
            for each of the six loci (e.g. ``{"A*1": "02:01", ...}``).
        :rtype: Dict[str, str]
        """
        # select two haplotypes based on the provided frequencies. `k=2`
        # means two draws, which can be the same (homozygous) or different
        chosen_haplotypes = random.choices(
            self.haplotypes, weights=self.frequencies, k=2
        )

        alleles1 = chosen_haplotypes[0].split('-')
        alleles2 = chosen_haplotypes[1].split('-')

        hla_type = {}
        for i, locus in enumerate(self.loci):
            hla_type[f"{locus}*1"] = alleles1[i]
            hla_type[f"{locus}*2"] = alleles2[i]

        return hla_type

    def generate_datasets(
            self,
            num_records: int = 100,
            match_ratio: float = 0.2,
            output_dir: str = '.',
            file_suffix: str = ''
    ) -> Dict[str, str]:
        """
        Generate and write recipient and donor CSV files to disk.

        Creates *num_records* recipients and *num_records* donors. A
        fraction of donors (controlled by *match_ratio*) are exact
        copies of randomly selected recipients, providing a known
        ground-truth for downstream validation.

        :param num_records: Number of records to generate for each of
            the recipient and donor files.
        :type num_records: int
        :param match_ratio: Fraction of donors that are exact HLA
            matches to a recipient, in the range ``[0.0, 1.0]``.
        :type match_ratio: float
        :param output_dir: Directory in which to write the output CSV
            files. Will be created if it does not exist.
        :type output_dir: str
        :param file_suffix: Suffix appended to output filenames before
            the ``.csv`` extension (e.g. ``"_1000"`` yields
            ``recipients_1000.csv``).
        :type file_suffix: str

        :returns: Ground-truth mapping of recipient IDs to their
            matching donor IDs (e.g.
            ``{"REC-0001": "DON-M-0001", ...}``).
        :rtype: Dict[str, str]
        """
        recipient_path = os.path.join(
            output_dir, f"recipients{file_suffix}.csv"
        )
        donor_path = os.path.join(
            output_dir, f"donors{file_suffix}.csv"
        )

        recipients = []
        donors = []
        match_ground_truth = {}

        # 1) generate all recipients
        for i in range(num_records):
            rec_id = f"REC-{i+1:04d}"
            hla_data = self._generate_individual()
            recipients.append({"ID": rec_id, **hla_data})

        # 2) generate donors, creating matches first
        num_matches = int(num_records * match_ratio)

        # create a shuffled list of recipients to pick from for matching
        recipients_to_match = recipients[:]
        random.shuffle(recipients_to_match)

        # create matching donors
        for i in range(num_matches):
            donor_id = f"DON-M-{i+1:04d}"
            matching_recipient = recipients_to_match[i]

            # create a donor with the exact same HLA data
            donor_hla = {
                key: val for key, val in matching_recipient.items()
                if key != "ID"
            }
            donors.append({"ID": donor_id, **donor_hla})

            # store the ground truth
            match_ground_truth[matching_recipient["ID"]] = donor_id

        # create non-matching donors
        num_non_matches = num_records - num_matches
        for i in range(num_non_matches):
            donor_id = f"DON-NM-{i+1:04d}"
            hla_data = self._generate_individual()
            donors.append({"ID": donor_id, **hla_data})

        # shuffle the donors list so matches are not all at the start
        random.shuffle(donors)

        # 3) write recipient csv
        rec_headers = ['ID'] + [
            f"patient-{locus}*{j}" for locus in self.loci for j in [1, 2]
        ]
        with open(recipient_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(rec_headers)
            for rec in recipients:
                row = [rec['ID']] + [
                    rec[f"{locus}*{j}"] for locus in self.loci
                    for j in [1, 2]
                ]
                writer.writerow(row)

        # 4) write donor csv
        don_headers = ['ID'] + [
            f"donor-{locus}*{j}" for locus in self.loci for j in [1, 2]
        ]
        with open(donor_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(don_headers)
            for don in donors:
                row = [don['ID']] + [
                    don[f"{locus}*{j}"] for locus in self.loci
                    for j in [1, 2]
                ]
                writer.writerow(row)

        print(f"Successfully generated '{recipient_path}' and '{donor_path}'.")
        print(f"Total Records: {num_records} per file.")
        print(f"Guaranteed Matches: {num_matches} ({match_ratio*100:.0f}%).")

        return match_ground_truth


if __name__ == "__main__":

    # haplotype frequency data
    haplotype_file_path = (
        "scripts/resources/dkms_haplotype_frequencies.txt"
    )
    with open(haplotype_file_path, 'r') as f:
        haplotype_data_string = f.read()

    generator = SyntheticDataGenerator(
        haplotype_freq_data=haplotype_data_string
    )

    dataset_sizes = [1000, 10000, 100000, 200000, 500000, 1000000]

    for num_records in dataset_sizes:
        ground_truth = generator.generate_datasets(
            num_records=num_records,
            match_ratio=0.9,
            output_dir="data/synthetic_data",
            file_suffix=f"_{num_records}"
        )
