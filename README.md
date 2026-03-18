# py-hla-match

## Regulatory Notice

The Py-HLA-Match open-source research software is not certified or conformity assessed as a medical device software or in-vitro medical device software and is intended to be used for research use only. It must therefore not be used for diagnosis or therapy of patients.

For more details on intended use, scope, and limitations, see the
[Software Card for Py-HLA-Match](SOFTWARE_CARD.md).

## About

The Py-HLA-Match open-source research software is a Python library for standardised, rule-based HLA (Human Leukocyte Antigen) matching in retrospective analyses, method development, benchmarking, and in-silico studies in immunogenetics and related fields.

## License

Copyright 2025 Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.

The Py-HLA-Match open-source research software is licensed under the Apache License, Version 2.0 (the "License").
You may obtain a copy of the License in the `LICENSE` file in this repository
or at http://www.apache.org/licenses/LICENSE-2.0.

## Terminology

The Py-HLA-Match open-source research software domain terms such as *patient*, *donor*, and *score* to mirror the structure of typical transplant research datasets (e.g. HSCT retrospective cohorts). These terms refer exclusively to roles and fields in research data and do **not** imply that Py-HLA-Match implements, recommends, or automates any clinical donor-selection or patient-management workflow.

All match levels and related outputs produced by the library are research metrics derived from HLA nomenclature semantics. They are **not** clinical risk scores or decision criteria.

![tests](https://github.com/tiadams/py-hla-match/actions/workflows/tests.yaml/badge.svg)

## Poetry Setup

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

```powershell
(Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | python -
```

## Installation

```bash
poetry install
```

## Testing

```bash
poetry run pytest
```
