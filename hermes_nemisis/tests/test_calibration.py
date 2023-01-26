import pytest
import hermes_nemisis.calibration as calib
from hermes_core.util.util import create_science_filename

level0_filename = "hermes_MAG_l0_2022339-000000_v0.bin"
level1_filename = "hermes_nms_l1_20221205_000000_v1.0.0.cdf"
ql_filename = "hermes_nms_ql_20221205_000000_v1.0.0.cdf"


def test_nemisis_sci_data_to_cdf():
    file_metadata = parse_science_filename(level0_filename)
    data = {}
    assert calib.nemisis_sci_data_to_cdf(data, file_metadata) == level1_filename


def test_calibrate_file():
    assert calib.calibrate_file(level0_filename) == level1_filename
    assert calib.calibrate_file(level1_filename) == ql_filename
    # with pytest.raises(ValueError) as excinfo:
    #    calib.calibrate_file("datafile_with_no_calib.cdf")
    # assert (
    #    str(excinfo.value)
    #    == "Calibration file for datafile_with_no_calib.cdf not found."
    # )


def test_process_file_level0():
    file_output = calib.process_file(level0_filename)
    assert len(file_output) == 1
    assert file_output[0] == level1_filename


def test_process_file_level1():
    file_output = calib.process_file(level1_filename)
    assert len(file_output) == 1
    assert file_output[0] == ql_filename


def test_get_calibration_file():
    assert calib.get_calibration_file("") is None


def test_read_calibration_file():
    assert calib.read_calibration_file("calib_file") is None
