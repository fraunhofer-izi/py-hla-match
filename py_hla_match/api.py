# api.py
"""
Simple API for HLA matching in allo-HSCT research.

Main function: py_hla_match()
"""

import logging
from pathlib import Path
from typing import Union, Dict, List, Optional, Tuple
import pandas as pd

from py_hla_match.loader import HLADataLoader
from py_hla_match.parser import HLAParser
from py_hla_match.matching import multi_locus_match
from py_hla_match.export import export_results
from py_hla_match.external import query_dpb1_tce

logger = logging.getLogger(__name__)


def py_hla_match(
    file_path: Union[str, Dict[str, str]],
    column_mapping: Dict[str, Dict[str, List[str]]],
    structure: str = 'paired',
    id_column: Optional[Union[str, Dict[str, str]]] = None,
    output_file: str = "hla_match_results.csv",
    dpb1_tce: bool = False
) -> pd.DataFrame:
    """
    Simple pipeline for HLA matching.

    Parameters
    ----------
    file_path : str or dict
        - str: Path file (for 'paired' or 'panel' structure)
        - dict: {'patient': 'path1.csv', 'donor': 'path2.csv'}
    column_mapping : dict
        Mapping from individuals to  HLA columns (by locus)
        Format: {
            'patient': {
                'A': ['col1', 'col2'],  # Two columns for A locus
                'B': ['col3', 'col4'],  # Two columns for B locus
                ...
            },
            'donor': {...}  # Same format
        }
    structure : {'paired', 'separate', 'panel'}
        - 'paired': Patient and donor in same file (default)
        - 'separate': Patient and donor in different files
        - 'panel': Multiple individuals
    id_column : str or dict, optional
        - str: Column name for IDs (for 'paired'/'panel')
        - dict: {'patient': 'id_col1', 'donor': 'id_col2'} (for 'separate')
        - None: Use row numbers (with warning)
    output_file : str
        Output filename
    dpb1_tce : bool
        Include DPB1 T-cell epitope matching analysis (default: False)

    Returns
    -------
    pd.DataFrame
        Matching results
        - patient_id, donor_id (or panel member IDs)
        - A_match_level1, A_match_level2, B_match_level1, ...
        - DPB1_TCE_status (if dpb1_tce=True)

    Examples
    --------
    # Paired data (most common)
    >>> results = py_hla_match(
    ...     file_path="transplant_data.xlsx",
    ...     column_mapping={
    ...         'patient': {
    ...             'A': ['Pat_A1', 'Pat_A2'],
    ...             'B': ['Pat_B1', 'Pat_B2'],
    ...             'C': ['Pat_C1', 'Pat_C2'],
    ...         },
    ...         'donor': {
    ...             'A': ['Don_A1', 'Don_A2'],
    ...             'B': ['Don_B1', 'Don_B2'],
    ...             'C': ['Don_C1', 'Don_C2'],
    ...         }
    ...     },
    ...     structure='paired',
    ...     id_column='PairID'
    ... )

    # Separate files
    >>> results = py_hla_match(
    ...     file_path={
    ...         'patient': 'patients.csv',
    ...         'donor': 'donors.csv'
    ...     },
    ...     column_mapping={...},
    ...     structure='separate',
    ...     id_column={
    ...         'patient': 'PatientID',
    ...         'donor': 'DonorID'
    ...     }
    ... )
    """
    # flatten
    flat_mapping = _flatten_column_mapping(column_mapping)

    # load based on structure
    if structure == 'paired':
        df = HLADataLoader.load(file_path, id_column=id_column)
        parser = HLAParser()
        pairs = parser.parse(df, 'paired', flat_mapping, id_column)
        results = _process_paired_matches(pairs, dpb1_tce)

    elif structure == 'separate':
        if not isinstance(file_path, dict) or 'patient' not in file_path \
                or 'donor' not in file_path:
            raise ValueError(
                "For 'separate' structure, file_path must be {'patient': "
                "'file1', 'donor': 'file2'}"
            )

        patient_id_col = id_column.get('patient') \
            if isinstance(id_column, dict) else None
        donor_id_col = id_column.get('donor') \
            if isinstance(id_column, dict) else None

        if not patient_id_col or not donor_id_col:
            logger.warning(
                "No ID columns specified for separate files. Using row numbers"
                " for matching."
            )

        patient_df = HLADataLoader.load(
            file_path['patient'], id_column=patient_id_col
        )
        donor_df = HLADataLoader.load(
            file_path['donor'], id_column=donor_id_col
        )

        # parse
        parser = HLAParser()
        patients = parser.parse(
            patient_df,
            'single',
            {'individual': flat_mapping['patient']},
            patient_id_col
        )
        donors = parser.parse(
            donor_df,
            'single',
            {'individual': flat_mapping['donor']},
            donor_id_col
        )

        # match by ID
        results = _process_separate_matches(patients, donors, dpb1_tce)

    elif structure == 'panel':
        df = HLADataLoader.load(file_path, id_column=id_column)
        parser = HLAParser()
        panels = parser.parse(df, 'panel', flat_mapping, id_column)
        results = _process_panel_matches(panels, dpb1_tce)

    else:
        raise ValueError(
            f"Unknown structure: {structure}. Must be 'paired', 'separate', "
            "or 'panel'"
        )

    # export
    export_results(
        results,
        output_file,
        errors=parser.errors if 'parser' in locals() else None,
        write_error_log=True
    )

    print(f"\nProcessed {len(results)} matches")
    print(f"Results saved to: {output_file}")

    error_file = Path(output_file).with_suffix('.errors.txt')
    if error_file.exists():
        print(f"Errors during parsing - see: {error_file}")

    return results


def _flatten_column_mapping(
        locus_mapping: Dict[str, Dict[str, List[str]]]
) -> Dict[str, List[str]]:
    """Convert locus-based mapping to flat list."""
    flat_mapping = {}

    for individual, loci in locus_mapping.items():
        columns = []
        for locus, locus_columns in loci.items():
            columns.extend(locus_columns)
        flat_mapping[individual] = columns

    return flat_mapping


def _process_paired_matches(
    pairs: List[Tuple[Union[str, int], object, object]],
    dpb1_tce: bool
) -> pd.DataFrame:
    """Process matching for paired data."""
    results_data = []
    dpb1_tce_query_required = False
    for pair_id, patient, donor in pairs:
        match_results = multi_locus_match(patient, donor)

        row_data = {
            'patient_id': pair_id,
            'donor_id': pair_id
        }

        # collect dpb1 for (potential) TCE
        dpb1_alleles = {}
        for match_result in match_results:
            locus = match_result.patient.locus
            if not locus:
                continue

            # Convert the match enum values to strings
            level1 = match_result.allele_match_levels[0].name
            level2 = match_result.allele_match_levels[1].name

            row_data[f'{locus}_match_level1'] = level1
            row_data[f'{locus}_match_level2'] = level2

            # If DPB1_TCE is enabled and this locus is DPB1...
            if dpb1_tce and locus == 'DPB1':
                # If both matches are ARD_MATCH or MATCH, skip TCE:
                if (
                    level1 in [
                        'ARD_MATCH',
                        'SYNONYMOUS_VARIANT_MATCH',
                        'NON_CODING_VARIANT_MATCH'
                    ]
                    and level2 in [
                        'ARD_MATCH',
                        'SYNONYMOUS_VARIANT_MATCH',
                        'NON_CODING_VARIANT_MATCH'
                    ]
                ):
                    # ARD match anyway
                    dpb1_tce_query_required = False
                    _extract_dpb1_alleles(match_result, dpb1_alleles)
                else:

                    dpb1_tce_query_required = True
                    _extract_dpb1_alleles(match_result, dpb1_alleles)

        if dpb1_tce and len(dpb1_alleles) == 4:
            if dpb1_tce_query_required:
                row_data['DPB1_TCE_status'] = \
                    _get_dpb1_tce_status(dpb1_alleles)
            else:
                row_data['DPB1_TCE_status'] = 'ARD_MATCH'

        results_data.append(row_data)

    return pd.DataFrame(results_data)


def _process_separate_matches(
    patients: List[Tuple[Union[str, int], object]],
    donors: List[Tuple[Union[str, int], object]],
    dpb1_tce: bool
) -> pd.DataFrame:
    """Process matching for separate files with shared IDs."""
    # ID lookup
    donor_dict = {donor_id: donor for donor_id, donor in donors}

    results_data = []
    for patient_id, patient in patients:
        if patient_id in donor_dict:
            donor = donor_dict[patient_id]
            match_results = multi_locus_match(patient, donor)

            row_data = {
                'patient_id': patient_id,
                'donor_id': patient_id
            }

            # copied from paired
            dpb1_alleles = {}
            for match_result in match_results:
                locus = match_result.patient.locus
                if locus:
                    row_data[f'{locus}_match_level1'] = \
                        match_result.allele_match_levels[0].name
                    row_data[f'{locus}_match_level2'] = \
                        match_result.allele_match_levels[1].name

                    if dpb1_tce and locus == 'DPB1':
                        _extract_dpb1_alleles(match_result, dpb1_alleles)

            if dpb1_tce and len(dpb1_alleles) == 4:
                row_data['DPB1_TCE_status'] = \
                    _get_dpb1_tce_status(dpb1_alleles)

            results_data.append(row_data)
        else:
            logger.warning(
                f"No matching donor found for patient ID: {patient_id}"
            )

    return pd.DataFrame(results_data)


def _process_panel_matches(
    panels: List[Tuple[Union[str, int], Dict[str, object]]],
    dpb1_tce: bool
) -> pd.DataFrame:
    """Process matching for panel data."""
    results_data = []

    for panel_id, members in panels:
        # first member is patient, others are potential donors
        member_names = list(members.keys())
        if len(member_names) < 2:
            logger.warning(f"Panel {panel_id} has only one member, skipping")
            continue

        patient_name = member_names[0]
        patient = members[patient_name]

        # match against members
        for donor_name in member_names[1:]:
            donor = members[donor_name]
            match_results = multi_locus_match(patient, donor)

            row_data = {
                'panel_id': panel_id,
                'patient': patient_name,
                'donor': donor_name
            }

            # Process match results
            dpb1_alleles = {}
            for match_result in match_results:
                locus = match_result.patient.locus
                if locus:
                    row_data[f'{locus}_match_level1'] = \
                        match_result.allele_match_levels[0].name
                    row_data[f'{locus}_match_level2'] = \
                        match_result.allele_match_levels[1].name

                    if dpb1_tce and locus == 'DPB1':
                        _extract_dpb1_alleles(match_result, dpb1_alleles)

            if dpb1_tce and len(dpb1_alleles) == 4:
                row_data['DPB1_TCE_status'] = \
                    _get_dpb1_tce_status(dpb1_alleles)

            results_data.append(row_data)

    return pd.DataFrame(results_data)


def _extract_dpb1_alleles(match_result, dpb1_alleles: dict) -> None:
    """Extract DPB1 alleles for TCE analysis."""
    if match_result.patient.hla1 and match_result.patient.hla1.allele:
        dpb1_alleles['patient_1'] = (
            f"{match_result.patient.hla1.allele_group}:"
            f"{match_result.patient.hla1.allele}"
        )
    if match_result.patient.hla2 and match_result.patient.hla2.allele:
        dpb1_alleles['patient_2'] = (
            f"{match_result.patient.hla2.allele_group}:"
            f"{match_result.patient.hla2.allele}"
        )
    if match_result.donor.hla1 and match_result.donor.hla1.allele:
        dpb1_alleles['donor_1'] = (
            f"{match_result.donor.hla1.allele_group}:"
            f"{match_result.donor.hla1.allele}"
        )
    if match_result.donor.hla2 and match_result.donor.hla2.allele:
        dpb1_alleles['donor_2'] = (
            f"{match_result.donor.hla2.allele_group}:"
            f"{match_result.donor.hla2.allele}"
        )


def _get_dpb1_tce_status(dpb1_alleles: dict) -> str:
    """Query DPB1 TCE API and return status."""
    try:
        required = ['patient_1', 'patient_2', 'donor_1', 'donor_2']
        for key in required:
            if key not in dpb1_alleles:
                logger.warning(f"Missing DPB1 allele: {key}")
                return "DPB1_DATA_INCOMPLETE"

        logger.debug(f"Querying DPB1 TCE with: {dpb1_alleles}")

        status = query_dpb1_tce(
            dpb1_alleles['patient_1'],
            dpb1_alleles['patient_2'],
            dpb1_alleles['donor_1'],
            dpb1_alleles['donor_2']
        )
        return status.value
    except Exception as e:
        logger.error(f"DPB1 TCE query failed: {e}")
        return "TCE_QUERY_FAILED"
