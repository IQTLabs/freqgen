from freqgen import generate

def test_1mer():
    assert generate({1: dict(A=0.5, T=0.5, G=0, C=0)}, "FK") == "TTTAAA"
