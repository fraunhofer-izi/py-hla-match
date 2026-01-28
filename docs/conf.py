# Configuration file for the Sphinx documentation builder.

from __future__ import annotations

import os
import sys


# Ensure the project root is on sys.path so autodoc can import the package.
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT_ROOT)


project = "py-hla-match"
author = "Tim Adams, Georg Popp"

# Keep version/release optional to avoid import side effects.
# If you prefer, we can read this from package metadata.

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autosummary_generate = True

autodoc_typehints = "description"

# If imports fail on Read the Docs due to optional deps, you can list them here.
autodoc_mock_imports: list[str] = []

templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"

