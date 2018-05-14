import pytest

from freqgen import gc_content

def test_gc_content():

    assert gc_content("ATGC") == 0.5
    assert gc_content("AAAAAATTTTTT") == 0.0
    assert gc_content("GCCC") == 1.0

def test_gc_content_base_checks():
    with pytest.raises(ValueError):
        gc_content("NOTAVALIDSEQUENCE")
