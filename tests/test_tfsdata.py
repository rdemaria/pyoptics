from pyoptics import tfsdata


def test_open():
    """Test opening a tfs table"""
    tfs = tfsdata.open("tests/tfsload.tfs")
    assert tfs["name"][0] == "NMCRING$START"
    assert tfs["s"][0] == 0.0
    assert tfs["name"][-1] == "NMCRING$END"
