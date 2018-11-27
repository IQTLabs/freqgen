import pytest

from freqgen import amino_acids_to_codons


def test_amino_acids_to_codons():
    assert amino_acids_to_codons("A", dict(GCA=1.0, GCT=0, GCG=0, GCC=0)) == "GCA"
    assert amino_acids_to_codons("AAA", dict(GCA=1.0, GCT=0, GCG=0, GCC=0)) == "GCAGCAGCA"

def test_bad_dict():
    with pytest.raises(ValueError):
        amino_acids_to_codons("AAA", dict(GCA=1.0, GCT=1.0, GCG=0, GCC=0))

    with pytest.raises(KeyError):
        amino_acids_to_codons("Q", dict(GCA=1.0, GCT=1.0, GCG=0, GCC=0))
