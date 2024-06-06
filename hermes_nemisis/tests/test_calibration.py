import pytest
from pathlib import Path

import hermes_nemisis.calibration as calib
from hermes_core.util.util import create_science_filename, parse_science_filename
from hermes_nemisis.calibration.readbin_CCSDS import readbin_CCSDS

hermes_nemisis_data_directory = Path('hermes_nemisis/data/')

level0_filename = "hermes_NEM_l0_2024085-174022_v00.bin"
level1_filename = "hermes_nem_l1_20240325T173118_v1.0.0.cdf"
ql_filename = "hermes_nem_ql_20240325T173118_v1.0.0.cdf"


@pytest.fixture(scope="session")
def level0_file(tmp_path_factory):
    fn = tmp_path_factory.mktemp("data") / level0_filename
    with open(fn, "w"):
        pass
    return fn


@pytest.fixture(scope="session")
def level1_file(tmp_path_factory):
    fn = tmp_path_factory.mktemp("data") / level1_filename
    with open(fn, "w"):
        pass
    return fn

# def test_readbin_CCSDS():  
#     output_file = readbin_CCSDS('hermes_nemisis/data/hermes_NEM_l0_2024085-174022_v00.bin')
#     # assert output_file.is_file()  # this fails.  output_file is a string object, not a file. is this a HermesData issue?

def test_parse_l0_sci_packets():
    data = calib.parse_l0_sci_packets('hermes_nemisis/data/hermes_NEM_l0_2024085-174022_v00.bin')
    assert data['MET_SECONDS'][0] == 2214
    assert len(data['MET_SECONDS']) == 270 

def test_nemisis_sci_data_to_l1(level0_file):
    """Test that the output filenames are correct and that a file was actually created."""
    data = calib.parse_l0_sci_packets('hermes_nemisis/data/hermes_NEM_l0_2024085-174022_v00.bin')
    output_file = calib.l0_sci_data_to_l1(data, level0_file)
    assert output_file.name == level1_filename
    assert output_file.is_file()


def test_calibrate_file_nofile_error():
    """Test that if file does not exist it produces the correct error. The file needs to be in the correct format."""
    with pytest.raises(FileNotFoundError):
        calib.calibrate_file(Path("hermes_NEM_l0_2032339-000000_v0.bin"))


def test_process_file_nofile_error():
    """Test that if file does not exist it produces the correct error. The file needs to be in the correct format."""
    with pytest.raises(FileNotFoundError):
        calib.process_file(Path("hermes_NEM_l0_2032339-000000_v0.bin"))


def test_calibrate_file(level0_file, level1_file):
    """Test that the output filenames are correct and that a file was actually created."""
#   NOTE the way the template is written, the simple file name is not sufficient
    output_file = calib.calibrate_file(level0_file)
    assert output_file.name == level1_filename
    assert output_file.is_file()
    output_file = calib.calibrate_file(level1_file)
    assert output_file.name == ql_filename
    assert output_file.is_file()

#     # with pytest.raises(ValueError) as excinfo:
#     #    calib.calibrate_file("datafile_with_no_calib.cdf")
#     # assert (
#     #    str(excinfo.value)
#     #    == "Calibration file for datafile_with_no_calib.cdf not found."
#     # )


# def test_process_file_level0(level0_file):
#     """Test that the output filenames are correct and that a file was actually created."""
#     file_output = calib.process_file(level0_file)
#     assert len(file_output) == 1
#     assert file_output[0].name == level1_filename
#     assert file_output[0].is_file()


def test_process_file_level1(level1_file):
    """Test that the output filenames are correct and that a file was actually created."""
    file_output = calib.process_file(level1_file)
    assert len(file_output) == 1
    assert file_output[0].name == ql_filename
    assert file_output[0].is_file()


def test_get_calibration_file():
    assert calib.get_calibration_file("") is None


def test_read_calibration_file():
    assert calib.read_calibration_file("calib_file") is None
