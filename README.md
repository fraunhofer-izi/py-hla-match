# py-hla-match

## Intended use and regulatory notice (non-license information)

Py-HLA-Match has been developed and evaluated **for research use**, such as
retrospective analyses, method development, benchmarking, and in‑silico studies in
immunogenetics and related fields.

Py-HLA-Match has **not** been developed, validated, or cleared for:

- diagnostic purposes
- clinical decision-making
- donor selection
- patient management
- real-time support of transplantation decisions

Py-HLA-Match is **not** a medical device or in vitro diagnostic medical device under
EU MDR/IVDR, FDA regulations, or any other medical device or IVD framework, and it has
not undergone regulatory approval or conformity assessment as such.

Any use of this software in a clinical, diagnostic, or patient‑related context is at
the sole responsibility of the user, who must ensure that all applicable legal and
regulatory requirements are identified, assessed, and fulfilled, including any
validation, risk management, and compliance obligations.

This section is provided for information and risk‑awareness only. 
**Nothing in this section limits, extends, or modifies the rights granted under the Apache License, Version 2.0, under which Py-HLA-Match is distributed.**


## Terminology

Py-HLA-Match uses domain terms such as *patient*, *donor* and *score* to
mirror the structure of typical transplant research datasets (e.g. HSCT
retrospective cohorts). These names indicate roles in research data and
do **not** imply that Py-HLA-Match implements, recommends, or automates
any clinical donor-selection or patient-management workflow.

All scores and match levels produced by the library are internal research
metrics derived from HLA nomenclature semantics. They are not validated or
intended as clinical risk scores or decision criteria.

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
