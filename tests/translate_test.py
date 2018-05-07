from freqgen import translate
import pytest

def test_standard_code():
    assert translate("ATTAATCAAACGGAGTTA") == "INQTEL"
    assert translate("GATACA") == "DT"

def test_bad_character():
    # invalid characters
    with pytest.raises(ValueError):
        translate("INQTEL")

    # invalid sequence length
    with pytest.raises(ValueError):
        translate("ATGGA")

def test_nonstandard_code():
    assert translate("CTC", genetic_code=3) == "T"
    assert translate("CTC") == "L"
