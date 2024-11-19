import pytest
from py_hla_match.main import hla_match_example

def test_hla_match_example():
    assert hla_match_example(1, 1) is True
    assert hla_match_example(1, 2) is False