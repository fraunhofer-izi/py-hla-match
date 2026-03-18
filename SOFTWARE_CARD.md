# Software Card for Py-HLA-Match

## Regulatory Notice

The Py-HLA-Match open-source research software is not certified or conformity 
assessed as a medical device software or in-vitro medical device software and 
is intended to be used for research use only. It must therefore not be used 
for diagnosis or therapy of patients. 

## Software Details

- **Name:** Py-HLA-Match (`py-hla-match`)
- **Version:** 0.1.0
- **Developed by:** Fraunhofer IZI, Fraunhofer SCAI
- **License:** Apache-2.0
- **Repository:** https://github.com/fraunhofer-izi/py-hla-match
- **PyPI package:** https://pypi.org/project/py-hla-match/

### Description

The Py-HLA-Match open-source research software is a Python library for standardized, rule-based computation of
HLA (Human Leukocyte Antigen) match categories. It reads HLA genotype data from
tabular formats, validates allele strings against HLA nomenclature, and
classifies pairwise allele comparisons into deterministic match/mismatch
categories.

No machine learning, statistical models, or trained parameters are used.
All results follow explicit, deterministic rules.

## Method Overview

1. Reads HLA genotype strings from tabular data (e.g. CSV, Excel)
2. Validates and parses allele strings against the HLA nomenclature into typed HLA objects
3. Assigns alleles to loci and pairs per individual
4. Computes pairwise match categories per locus (e.g. Allele Match, Antigen Mismatch, ARD Match) using deterministic, rule-based logic
5. Optional: performs DPB1 T-cell epitope classification via the EBI web service
6. Returns tabular output (e.g. pandas DataFrame / CSV) with match categories per locus

## Risks and Limitations

Match categories depend on the configured IMGT/HLA database version. Users must ensure that the chosen version iscorrect  for their research data.

The optional DPB1 lookup depends on an external service (EBI). Availability and behaviour of this service are outside the control of this library

Users are responsible for assessing and validating the suitability of the Py-HLA-Match open-source research software in their specific context.

## Evaluation

Correctness is tested with a deterministic test suite covering:

- allele parsing and validation
- core allele matching logic
- allele pairing behaviour
- configuration
- edge cases

## Technical Profile

- **Language:** Python ≥ 3.9
- **Runtime dependencies:** pandas, py-ard, openpyxl, requests

## Contact

- Tim Adams: tim.adams@scai.fraunhofer.de (Fraunhofer SCAI)
- Georg Popp: georg.popp@izi.fraunhofer.de (Fraunhofer IZI)