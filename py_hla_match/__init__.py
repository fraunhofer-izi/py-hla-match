# __init__.py
"""
py-hla-match: A package for HLA matching.
"""

__version__ = "0.1.0"

from .hla import HLA

from .exceptions import MalformedHLAStringError

__all__ = [
    'hla',
    'exceptions',
    'matching'

]