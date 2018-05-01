import pytest
from seqgen import amino_acid_seq

def test_length():
    assert len(amino_acid_seq(length=10, frequencies=dict(A=1))) == 10
    with pytest.raises(ValueError):
        amino_acid_seq(length=0, frequencies=dict(A=1))
    with pytest.raises(ValueError):
        amino_acid_seq(length=-1, frequencies=dict(A=1))

def test_frequencies_sum():
    with pytest.raises(ValueError):
        amino_acid_seq(length=1, frequencies=dict(A=0.5))
    assert amino_acid_seq(length=5, frequencies=dict(A=1)) == "AAAAA"
