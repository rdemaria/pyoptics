from pyoptics import optics


def test_optics_open():
    """Test opening an optics file"""
    t = optics.open("tests/tfsload.tfs")
    assert t.name[0] == "NMCRING$START"
    assert t.s[0] == 0.0
    assert t.name[-1] == "NMCRING$END"


def test_optics_header():
    """Test reading the header of an optics file"""
    t = optics.open("tests/tfsload.tfs")
    assert t.header["name"] == "TWISS"
