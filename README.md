# Py-HLA-Match

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19593513.svg)](https://doi.org/10.5281/zenodo.19593513)
![tests](https://github.com/fraunhofer-izi/py-hla-match/actions/workflows/tests.yaml/badge.svg)
[![codecov](https://codecov.io/gh/fraunhofer-izi/py-hla-match/branch/main/graph/badge.svg)](https://codecov.io/gh/fraunhofer-izi/py-hla-match)
![docs](https://github.com/fraunhofer-izi/py-hla-match/actions/workflows/docs.yaml/badge.svg)
![version](https://img.shields.io/pypi/v/py-hla-match)
![license](https://img.shields.io/badge/license-Apache%202.0-blue)

## About

Py-HLA-Match is a Python library for standardised, rule-based HLA (Human Leukocyte Antigen) matching in retrospective analyses, method development, benchmarking, and in-silico studies in immunogenetics and related fields.

## Regulatory Notice

Py-HLA-Match is **not** certified or conformity assessed as a medical device software or in-vitro medical device software and is intended for **research use only**. It must therefore not be used for diagnosis or therapy of patients.

For more details on intended use, scope, and limitations, see the [Software Card](SOFTWARE_CARD.md).

## Installation

Install from PyPI:

```bash
pip install py-hla-match
```

## Quickstart

This quickstart uses the artificial CSVs bundled under the `demo` folder and avoids any real or sensitive data.

### Run a basic pairwise match

Use the synthetic patient and donor CSVs and write results to a new file:

```python
from py_hla_match.parser import HLADataSource
from py_hla_match.export import PairwiseMatch

data_path = "py_hla_match/demo/data/random_data/synthetic_patients.csv"
donor_path = "py_hla_match/demo/data/random_data/synthetic_donors.csv"
output_path = "py_hla_match/demo/data/random_data/match_results.csv"

src = HLADataSource(
    data_path,
    col_idx_start=1,
    col_idx_stop=13,
    row_idx_start=1,
)

tgt = HLADataSource(
    donor_path,
    col_idx_start=1,
    col_idx_stop=13,
    row_idx_start=1,
)

matcher = PairwiseMatch(
    source=src,
    target=tgt,
    storage_filename=output_path,
    resolution="high",
)

matcher.run()
```

### Inspect raw allele-level results

Convert raw match levels to a DataFrame and write to CSV:

```python
raw_output_path = "py_hla_match/demo/data/random_data/match_results_raw.csv"
matcher.raw_to_df().to_csv(raw_output_path, index=False)
```

## Terminology

Py-HLA-Match uses domain terms such as *patient*, *donor*, and *score* to mirror the structure of typical transplant research datasets (e.g. HSCT retrospective cohorts). These terms refer exclusively to roles and fields in research data and do **not** imply that Py-HLA-Match implements, recommends, or automates any clinical donor-selection or patient-management workflow.

All match levels and related outputs produced by the library are research metrics derived from HLA nomenclature semantics. They are **not** clinical risk scores or decision criteria.

## Development

### Prerequisites

Install [Poetry](https://python-poetry.org/docs/#installation):

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

On Windows (PowerShell):

```powershell
(Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | python -
```

### Setup

```bash
git clone https://github.com/fraunhofer-izi/py-hla-match.git
cd py-hla-match
poetry install
```

### Running tests

```bash
poetry run pytest
```

## License

Copyright 2025 Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.

Licensed under the Apache License, Version 2.0. You may obtain a copy of the License in the [`LICENSE`](LICENSE) file or at <http://www.apache.org/licenses/LICENSE-2.0>.
