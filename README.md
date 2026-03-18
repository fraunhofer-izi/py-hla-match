# py-hla-match

## Intended use and regulatory notice (non-license information)

Py-HLA-Match is research software intended **for research use**, such as retrospective analyses, method development, benchmarking, and in‑silico studies in immunogenetics and related fields. Py-HLA-Match is intended to be used on research data.

Py-HLA-Match is **not** a medical device or in vitro diagnostic medical device under EU MDR/IVDR, FDA regulations, or any comparable medical device or IVD regulations, and it has **not** undergone regulatory approval or conformity assessment as such. No conformity assessment, CE marking, 510(k), De Novo, PMA, or other regulatory approval has been sought or obtained. Py-HLA-Match is **not** intended for clinical diagnosis, donor selection, patient management, or for supporting or influencing medical or transplantation decisions for individual patients.

Any use of this software in a clinical, diagnostic, or patient‑related context is at the sole responsibility of the user or deploying organisation, including compliance with any applicable legal and regulatory requirements.

This section is provided solely for information and risk awareness. Nothing in this section limits, extends, or modifies the rights granted under the Apache License, Version 2.0, under which Py‑HLA‑Match is distributed.


## Terminology

Py-HLA-Match uses domain terms such as *patient*, *donor* and *_score* to mirror the structure of typical transplant research datasets (e.g. HSCT
retrospective cohorts). These terms refer exclusively to roles and fields in research data and do **not** imply that Py-HLA-Match implements, recommends, or automates any clinical donor-selection or patient-management workflow.

All scores, match levels, and related outputs produced by the library are intended as research metrics derived from HLA nomenclature semantics. They are **not** intended as clinical risk scores or decision criteria.

![tests](https://github.com/tiadams/py-hla-match/actions/workflows/tests.yaml/badge.svg)
![docs](https://github.com/tiadams/py-hla-match/actions/workflows/docs.yaml/badge.svg)
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
