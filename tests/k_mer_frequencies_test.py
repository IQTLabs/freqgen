from freqgen import k_mer_frequencies, codon_frequencies
import numpy as np
import pytest
from hypothesis import assume, given, strategies as st
from itertools import product

def test_amino_acid():
    assert k_mer_frequencies("INQTEL", 1, include_missing=False) == {1: {'E': 0.16666666666666666,
                                                                         'I': 0.16666666666666666,
                                                                         'L': 0.16666666666666666,
                                                                         'N': 0.16666666666666666,
                                                                         'Q': 0.16666666666666666,
                                                                         'T': 0.16666666666666666}}

def test_dna():
    assert k_mer_frequencies("GATGATGGC", 3, include_missing=False) == {3: {'ATG': 0.2857142857142857,
                                                                            'GAT': 0.2857142857142857,
                                                                            'GGC': 0.14285714285714285,
                                                                            'TGA': 0.14285714285714285,
                                                                            'TGG': 0.14285714285714285}}

def test_include_missing():
    assert k_mer_frequencies("GATGATGGC", 2, include_missing=True) == {2: {'AA': 0,
                                                                           'AC': 0,
                                                                           'AG': 0,
                                                                           'AT': 0.25,
                                                                           'CA': 0,
                                                                           'CC': 0,
                                                                           'CG': 0,
                                                                           'CT': 0,
                                                                           'GA': 0.25,
                                                                           'GC': 0.125,
                                                                           'GG': 0.125,
                                                                           'GT': 0,
                                                                           'TA': 0,
                                                                           'TC': 0,
                                                                           'TG': 0.25,
                                                                           'TT': 0}}
@given(st.text(alphabet="ATGCatgc"))
def test_k_values(s):
    assume(len(s) > 6)
    # make sure that the lengths are right
    for i in range(1, 6):
        assert len(k_mer_frequencies(s, i, include_missing=True)[i]) == 4**i
        assert len(k_mer_frequencies(s, i, include_missing=True, vector=True)) == 4**i

    # ensure invalid values of k raise an error
    with pytest.raises(ValueError):
        k_mer_frequencies(s, 0)
    with pytest.raises(ValueError):
        k_mer_frequencies(s, -1)

def test_vectorization():
    # check that the ordering is alphabetical
    assert np.array_equal(k_mer_frequencies("A", 1, include_missing=True, vector=True), np.array([1.0, 0.0, 0.0, 0.0]))
    assert np.array_equal(k_mer_frequencies("T", 1, include_missing=True, vector=True), np.array([0.0, 0.0, 0.0, 1.0]))
    assert np.array_equal(k_mer_frequencies("G", 1, include_missing=True, vector=True), np.array([0.0, 0.0, 1.0, 0.0]))
    assert np.array_equal(k_mer_frequencies("C", 1, include_missing=True, vector=True), np.array([0.0, 1.0, 0.0, 0.0]))
    assert np.array_equal(k_mer_frequencies("AT", 1, include_missing=True, vector=True), np.array([0.5, 0.0, 0.0, 0.5]))

def test_multiple_k():
    assert np.array_equal(k_mer_frequencies("AA", [1, 2], include_missing=True, vector=True), np.array([1.0, 0.0, 0.0, 0.0, # k = 1
                                                                                                        1.0, 0.0, 0.0, 0.0, # k = 2
                                                                                                        0.0, 0.0, 0.0, 0.0,
                                                                                                        0.0,  0.0, 0.0, 0.0,
                                                                                                        0.0,  0.0, 0.0, 0.0]))

    assert np.array_equal(k_mer_frequencies("AA", [2, 1], include_missing=True, vector=True), np.array([1.0, 0.0, 0.0, 0.0, # k = 1
                                                                                                        1.0, 0.0, 0.0, 0.0, # k = 2
                                                                                                        0.0, 0.0, 0.0, 0.0,
                                                                                                        0.0,  0.0, 0.0, 0.0,
                                                                                                        0.0,  0.0, 0.0, 0.0]))
    assert k_mer_frequencies("AA", [1, 2], include_missing=False) == {1: {"A": 1.0},
                                                                      2: {"AA": 1.0}}

def test_multiple_seqs():
    assert k_mer_frequencies(["A", "A"], 1, include_missing=False) == {1: {"A": 1.0}}
    assert k_mer_frequencies(["A", "T"], 1, include_missing=False) == {1: {"A": 0.5, "T": 0.5}}
    assert k_mer_frequencies(["AA", "TT"], 2, include_missing=False) == {2: {"AA": 0.5, "TT": 0.5}}
    assert k_mer_frequencies(["A", "T"], 1, include_missing=True) == {1: {"A": 0.5, "T": 0.5, "G": 0.0, "C": 0.0}}
    assert np.array_equal(k_mer_frequencies(["A", "T"], 1, include_missing=True, vector=True), np.array([0.5, 0.0, 0.0, 0.5]))

def test_invalid_args():
    with pytest.raises(ValueError):
        k_mer_frequencies("A", 1, include_missing=False, vector=True)
    with pytest.raises(ValueError):
        k_mer_frequencies("", 1)

@given(st.lists(elements=st.sampled_from(["".join(codon) for codon in product("ATGC", repeat=3)]), min_size=1))
def test_codon_frequencies(s):
    assert k_mer_frequencies(s, 1, codons=True)["codons"] == codon_frequencies(s)
