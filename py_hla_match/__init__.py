# __init__.py
"""
py-hla-match: A package for HLA matching.
"""

__version__ = "0.1.0"

__all__ = [
    'hla',
    'exceptions',
    'matching',
    'models',
    'entities'
]

# import global ARD instance
from .singleton import get_ard_instance

global_ard = get_ard_instance()
