# Configuration file for the Sphinx documentation builder.

from __future__ import annotations

import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath("../.."))

project = "py-hla-match"
author = "Georg Popp, Tim Adams"
copyright = f"{datetime.now().year}, {author}"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autosummary_generate = True
autosummary_generate_overwrite = True

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Keep module descriptions in docstrings for readability.
autodoc_typehints = "description"
