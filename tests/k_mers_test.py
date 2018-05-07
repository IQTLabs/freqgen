from freqgen import k_mers
from types import GeneratorType

def test_k_mers():
    assert list(k_mers("GATTACA", 1)) == ['G', 'A', 'T', 'T', 'A', 'C', 'A']
    assert list(k_mers("GATTACA", 2)) == ['GA', 'AT', 'TT', 'TA', 'AC', 'CA']
    assert list(k_mers("GATTACA", 3)) == ['GAT', 'ATT', 'TTA', 'TAC', 'ACA']
    assert list(k_mers("GATTACA", 4)) == ['GATT', 'ATTA', 'TTAC', 'TACA']
    assert type(k_mers("GATTACA", 4)) == GeneratorType
