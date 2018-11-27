from types import GeneratorType

import pytest

from freqgen import k_mers


def test_k_mers():
    assert list(k_mers("GATTACA", 1)) == ['G', 'A', 'T', 'T', 'A', 'C', 'A']
    assert list(k_mers("GATTACA", 2)) == ['GA', 'AT', 'TT', 'TA', 'AC', 'CA']
    assert list(k_mers("GATTACA", 3)) == ['GAT', 'ATT', 'TTA', 'TAC', 'ACA']
    assert list(k_mers("GATTACA", 4)) == ['GATT', 'ATTA', 'TTAC', 'TACA']
    assert type(k_mers("GATTACA", 4)) == GeneratorType

def test_length_check():
    with pytest.raises(ValueError):
        list(k_mers("A", 4))
    with pytest.raises(ValueError):
        list(k_mers("", 0))
    with pytest.raises(ValueError):
        list(k_mers("A", 0))
    with pytest.raises(ValueError):
        list(k_mers("", 4))
