"""
A module for all things calibration.
"""
import random
import os.path

from spacepy import pycdf

import ccsdspy
from ccsdspy import PacketField

from hermes_core import log
from hermes_core.util.util import create_science_filename, parse_science_filename

import hermes_nemisis
from hermes_nemisis.io import read_file

__all__ = ["calibrate_file", "get_calibration_file", "read_calibration_file"]


def process_file(data_filename: str) -> list:
    """
    This is the entry point for the pipeline processing.
    It runs all of the various processing steps required.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of an input file

    Returns
    -------
    output_filenames: list
        Fully specificied filenames for the output files.
    """
    log.info(f"Processing file {data_filename}.")

    output_files = []

    calibrated_file = calibrate_file(data_filename)
    output_files.append(calibrated_file)
    #  data_plot_files = plot_file(data_filename)
    #  calib_plot_files = plot_file(calibrated_file)

    # add other tasks below
    return output_files


def calibrate_file(data_filename) -> str:
    """
    Given an input data file, raise it to the next level
    (e.g. level 0 to level 1, level 1 to quicklook) it and return a new file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of the input data file.

    Returns
    -------
    output_filename: str
        Fully specificied filename of the output file.

    Examples
    --------
    >>> from hermes_nemisis.calibration import calibrate_file
    >>> level1_file = calibrate_file('hermes_MAG_l0_2022239-000000_v0.bin')  # doctest: +SKIP
    """
    log.info(f"Calibrating file:{data_filename}.")
    output_filename = (
        data_filename  # for testing, the output filename MUST NOT same as input
    )
    file_metadata = parse_science_filename(data_filename)

    # check if level 0 binary file, if so call appropriate functions
    if (
        file_metadata["instrument"] == hermes_nemisis.INST_NAME
        and file_metadata["level"] == "l0"
    ):
        data = parse_nemisis_sci_packets(data_filename)
        output_filename = nemisis_sci_data_to_cdf(data, file_metadata)
    elif file_metadata["instrument"] == "sunsensor" and file_metadata["level"] == "l1":
        data = read_file(data_filename)
        calib_file = get_calibration_file(data_filename)
        if calib_file is None:
            raise ValueError(f"Calibration file for {data_filename} not found.")
        else:
            calib_data = read_calibration_file(calib_file)
        # calibrate your data
        # create a new file name with create_science_filename()
        # write a new cdf file

        # example log messages
        log.info(f"Despiking removing {random.randint(0, 10)} spikes")
        log.warning(f"Despiking could not remove {random.randint(1, 5)}")

    else:
        raise ValueError(f"The file {data_filename} is not recognized.")

    return output_filename


def parse_nemisis_sci_packets(data_filename) -> dict:
    """
    Parse a level 0 nemisis binary file containing CCSDS packets.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename

    Returns
    -------
    result: dict
        A dictionary of arrays which includes the ccsds header fields

    Examples
    --------
    >>> import hermes_nemisis.calibration as calib
    >>> data_filename = "hermes_MAG_l0_2022339-000000_v0.bin"
    >>> data = calib.parse_nemisis_sci_packets(data_filename)  # doctest: +SKIP
    """
    log.info(f"Parsing packets from file:{data_filename}.")

    pkt = ccsdspy.FixedLength.from_file(
        os.path.join(hermes_nemisis._data_directory, "MAG_sci_packet_def.csv")
    )
    data = pkt.load(data_filename)
    return data


def nemisis_sci_data_to_cdf(data, file_metadata) -> str:
    """
    Write level 0 nemisis science data to a level 1 cdf file.

    Parameters
    ----------
    data: dict
        A dictionary of arrays which includes the ccsds header fields
    metadata: dict
        A metadata dictionary from the originating binary file.
        The output from `parse_science_filename()`.

    Returns
    -------
    output_filename: str
        Fully specificied filename of cdf file

    Examples
    --------
    >>> from hermes_core.util.util import parse_science_filename
    >>> import hermes_nemisis.calibration as calib
    >>> data_filename = "hermes_MAG_l0_2022339-000000_v0.bin"
    >>> metadata = parse_science_filename(data_filename)  # doctest: +SKIP
    >>> data_packets = calib.parse_nemisis_sci_packets(data_filename)  # doctest: +SKIP
    >>> cdf_filename = calib.nemisis_sci_data_to_cdf(data_packets, metadata)  # doctest: +SKIP
    """

    cdf_filename = create_science_filename(
        hermes_nemisis.INST_NAME, file_metadata["time"], "l1", "1.0.0", test=False
    )

    return cdf_filename


def get_calibration_file(data_filename, time=None) -> str:
    """
    Given a time, return the appropriate calibration file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)
    time: ~astropy.time.Time

    Returns
    -------
    calib_filename: str
        Fully specificied filename for the appropriate calibration file.

    Examples
    --------
    """
    return None


def read_calibration_file(calib_filename):
    """
    Given a calibration, return the calibration structure.

    Parameters
    ----------
    calib_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)

    Returns
    -------
    output_filename: str
        Fully specificied filename of the appropriate calibration file.

    Examples
    --------
    """

    # if can't read the file

    return None
