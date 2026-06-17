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

### Run hla matching (pairwise := )

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
    include_ard_details=True,
    include_molecular_details=True,
    include_dpb1_tce=False,
    include_homozygosity=False,
    overwrite=True,
)

matcher.run()
```

### Inspect raw allele-level results

Convert raw match levels to a DataFrame and write to CSV:

```python
df = matcher.to_df()
print(df.head())
```

## Matching Logic

Py-HLA-Match classifies each donor–recipient allele pair through in two stages
that follow IPD-IMGT/HLA nomenclature semantics.

### Stage 1: Mismatch Detection

Both alleles are reduced to their ARD (antigen recognition domain) equivalent
via P-group affiliation. If ARD representations differ:

| Condition | Classification |
| ----------- | --------------- |
| Field 1 (allele group) differs | `ANTIGEN_MISMATCH` |
| Field 1 identical, ARD field 2 differs | `ALLELE_MISMATCH` |
| DRB3/4/5: same broad locus, different sublocus | `SUBLOCUS_MISMATCH` |
| Insufficient resolution for comparison | `NOT_ASSESSABLE` |

### Stage 2: Match Refinement (ARD-matched pairs only)

Pairs classified as `ARD_MATCH` are refined along two independent dimensions:

**ARD match level**: identity at the antigen recognition domain:

| Level | Meaning |
| ------- | --------- |
| `P_GROUP_MATCH` | Identical ARD amino acid sequence |
| `G_GROUP_MATCH` | Identical ARD nucleotide sequence |

**Molecular match level**: depth of sequence identity beyond ARD:

| Level | Condition | Example |
| ------- | ----------- | --------- |
| `NOT_ASSESSABLE` | Only 2 fields typed | `A*02:01` vs `A*02:01` |
| `ARD_MATCH_ONLY` | Field 2 fields differ but allels share P group | `A*01:01` vs `A*01:510` |
| `FULL_PROTEIN_MATCH` | Fields 1–2 identical, field 3 differs or untyped | `A*02:01:01` vs `A*02:01:02` |
| `CODING_SEQUENCE_MATCH` | Fields 1–3 identical, field 4 differs or untyped | `A*02:01:01:01` vs `A*02:01:01:02` |
| `EXACT_ALLELE_MATCH` | All 4 fields identical | `A*02:01:01:01` vs `A*02:01:01:01` |

### Certainty

Each level carries a certainty indicator:

- **`CERTAIN`**: typing resolution is sufficient to confirm the level
- **`UNCERTAIN`**: a higher level remains possible given untyped fields

## Examples

The examples below illustrate key design choices that Py-HLA-Match makes explicit.

Each is drawn directly from the test suite and is independently reproducible.

### Resolution-aware certainty

The same ARD match is classified differently depending on typing depth:

```python
from py_hla_match.hla import HLA
from py_hla_match.matching import allele_pair_match
from py_hla_match.models import HLAPair

# 4-field identical -> EXACT_ALLELE_MATCH, CERTAIN
patient = HLAPair(HLA("A*01:01:01:01"), HLA("A*02:01:01:01"))
donor   = HLAPair(HLA("A*02:01:01:01"), HLA("A*01:01:01:01"))
r = allele_pair_match(patient, donor)
# r.molecular_match_levels      -> (EXACT_ALLELE_MATCH, EXACT_ALLELE_MATCH)
# r.molecular_match_certainties -> (CERTAIN, CERTAIN)

# 4-field, field 4 differs -> CODING_SEQUENCE_MATCH, CERTAIN
patient = HLAPair(HLA("A*01:01:01:01"), HLA("A*01:01:01:04"))
donor   = HLAPair(HLA("A*01:01:01:03"), HLA("A*01:01:01:05"))
r = allele_pair_match(patient, donor)
# r.molecular_match_levels      -> (CODING_SEQUENCE_MATCH, CODING_SEQUENCE_MATCH)
# r.molecular_match_certainties -> (CERTAIN, CERTAIN)

# 3-field vs 4-field -> CODING_SEQUENCE_MATCH, UNCERTAIN
patient = HLAPair(HLA("A*01:01:01"), HLA("A*01:02:01"))
donor   = HLAPair(HLA("A*01:02:01:01"), HLA("A*01:01:01:03"))
r = allele_pair_match(patient, donor)
# r.molecular_match_levels      -> (CODING_SEQUENCE_MATCH, CODING_SEQUENCE_MATCH)
# r.molecular_match_certainties -> (UNCERTAIN, UNCERTAIN)

# 2-field identical -> FULL_PROTEIN_MATCH, UNCERTAIN
patient = HLAPair(HLA("A*01:01"), HLA("A*01:02"))
donor   = HLAPair(HLA("A*01:02"), HLA("A*01:01"))
r = allele_pair_match(patient, donor)
# r.molecular_match_levels      -> (FULL_PROTEIN_MATCH, FULL_PROTEIN_MATCH)
# r.molecular_match_certainties -> (UNCERTAIN, UNCERTAIN)
```

### ARD equivalence is not sequence identity

Alleles with different names can share the same ARD reduction. Py-HLA-Match
explicitly distinguishes immunological equivalence from sequence identity:

```python
# A*02:01 and A*02:09 share the same G-group but differ at field 2
patient = HLAPair(HLA("A*02:01:01"), HLA("A*02:01:01"))
donor   = HLAPair(HLA("A*02:09:01"), HLA("A*02:09:01"))
r = allele_pair_match(patient, donor)
# r.allele_match_levels    -> (ARD_MATCH, ARD_MATCH)
# r.ard_match_levels       -> (G_GROUP_MATCH, G_GROUP_MATCH)
# r.molecular_match_levels -> (ARD_MATCH_ONLY, ARD_MATCH_ONLY)
```

`ARD_MATCH_ONLY` indicates that the alleles are equivalent at the antigen
recognition domain but differ in their full sequence.

### Expression suffix policy

Expression suffixes (N, L, S, C, A, Q) are evaluated under a configurable
policy. The default treats risk-associated suffixes as functional mismatches
and questionable expression (Q) as not assessable:

```python
from py_hla_match.matching import allele_match

# Null allele vs expressed -> ALLELE_MISMATCH (default)
allele_match(HLA("C*03:693"), HLA("C*03:20N"))   # -> ALLELE_MISMATCH

# Questionable expression -> NOT_ASSESSABLE (default)
allele_match(HLA("A*01:436Q"), HLA("A*01:01:70"))  # -> NOT_ASSESSABLE
```

Override the default to apply center-specific conventions:

```python
from py_hla_match.policy import ExpressionSuffixPolicy, ExpressionSuffixMatchLevel
from py_hla_match.config import HLAMatchConfig, set_config

set_config(HLAMatchConfig(expression_suffix_policy=ExpressionSuffixPolicy(
    q_present=ExpressionSuffixMatchLevel.ALLELE_MISMATCH,
)))
allele_match(HLA("A*01:436Q"), HLA("A*01:01:70"))  # -> ALLELE_MISMATCH
```

### DRB3/4/5 sublocus mismatch

The DRB3/4/5 region is normalized to a shared DRB345 locus and is given an additional mismatch
class:

```python
# Different subloci within DRB345
allele_match(HLA("DRB3*02:02:01"), HLA("DRB4*01:03:01"))  # -> DRB345_SUBLOCUS_MISMATCH

# Present sublocus vs non-expressed marker
allele_match(HLA("DRB3*01"), HLA("DRBX*NE"))              # -> DRB345_SUBLOCUS_MISMATCH
```

### Insufficient resolution

When typing resolution is too low for ARD comparison, the result is explicitly
flagged rather than silently excluded or assumed to match:

```python
# 1-field cannot confirm ARD equivalence even within the same allele group
allele_match(HLA("B*07"), HLA("B*07:05"))  # -> NOT_ASSESSABLE

# Missing data
allele_match(HLA("A*NE"), HLA("A*01:01"))  # -> NOT_ASSESSABLE
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
