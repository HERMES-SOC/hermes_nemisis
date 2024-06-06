"""
A module for all things calibration.
"""

import random
from pathlib import Path

import ccsdspy

from hermes_core import log
from hermes_core.util.util import create_science_filename, parse_science_filename

import hermes_nemisis
from hermes_nemisis.io import read_file

__all__ = [
    "process_file",
    "parse_l0_sci_packets",
    "l0_sci_data_to_l1",
    "calibrate_file",
    "get_calibration_file",
    "read_calibration_file",
]


def process_file(data_filename: Path) -> list:
    """
    This is the entry point for the pipeline processing.
    It runs all of the various processing steps required.

    Parameters
    ----------
    data_filename: `pathlib.Path`
        Fully specificied filename of an input file.
        The file contents: Traditional binary packets in CCSDS format

    Returns
    -------
    output_filenames: `list[pathlib.Path]`
        Fully specificied filenames for the output CDF files.
    """
    log.info(f"Processing file {data_filename}.")
    output_files = []

    calibrated_file = calibrate_file(data_filename)
    output_files.append(calibrated_file)
    #  data_plot_files = plot_file(data_filename)
    #  calib_plot_files = plot_file(calibrated_file)

    # add other tasks below
    return output_files


def calibrate_file(data_filename: Path) -> Path:
    """
    Given an input data file, raise it to the next level
    (e.g. level 0 to level 1, level 1 to quicklook) it and return a new file.

    Parameters
    ----------
    data_filename: `pathlib.Path`
        Fully specificied filename of the input data file.

    Returns
    -------
    output_filename: `pathlib.Path`
        Fully specificied filename of the output file.

    Examples
    --------
    >>> from hermes_nemisis.calibration import calibrate_file
    >>> level1_file = calibrate_file('hermes_NEM_l0_2022239-000000_v0.bin')  # doctest: +SKIP
    """
    log.info(f"Calibrating file:{data_filename}.")
    file_metadata = parse_science_filename(data_filename.name)

    # check if level 0 binary file, if so call appropriate functions
    if (
        file_metadata["instrument"] == hermes_nemisis.INST_NAME
        and file_metadata["level"] == "l0"
    ):
        data = parse_l0_sci_packets(data_filename)
        level1_filename = l0_sci_data_to_l1(data, data_filename)
        output_filename = level1_filename
    elif (
        file_metadata["instrument"] == hermes_nemisis.INST_NAME
        and file_metadata["level"] == "l1"
    ):
        # generate the quicklook data
        #
        # the following shows an example flow for calibrating a file
        # data = read_file(data_filename)
        # calib_file = get_calibration_file(data_filename)
        # if calib_file is None:
        #    raise ValueError(f"Calibration file for {data_filename} not found.")
        # else:
        #    calib_data = read_calibration_file(calib_file)

        # test opening the file
        with open(data_filename, "r") as fp:
            pass

        # now that you have your calibration data, you can calibrate the science data
        ql_filename = data_filename.parent / create_science_filename(
            file_metadata["instrument"],
            file_metadata["time"],
            "ql",
            file_metadata["version"],
        )

        # write your cdf file below
        # create an empty file for testing purposes
        with open(data_filename.parent / ql_filename, "w"):
            pass

        # example log messages
        log.info(f"Despiking removing {random.randint(0, 10)} spikes")
        log.warning(f"Despiking could not remove {random.randint(1, 5)}")
        output_filename = ql_filename
    else:
        raise ValueError(f"The file {data_filename} is not recognized.")

    return output_filename


def parse_l0_sci_packets(data_filename: Path) -> dict:
    """
    Parse a level 0 nemisis binary file containing CCSDS packets.

    Parameters
    ----------
    data_filename: `pathlib.Path`
        Fully specificied filename

    Returns
    -------
    result: dict
        A dictionary of arrays which includes the ccsds header fields

    Examples
    --------
    >>> import hermes_nemisis.calibration as calib
    >>> data_filename = "hermes_NEM_l0_2022339-000000_v0.bin"
    >>> data = calib.parse_nemisis_sci_packets(data_filename)  # doctest: +SKIP
    """
    log.info(f"Parsing packets from file:{data_filename}.")

    pkt = ccsdspy.FixedLength.from_file(
        Path(hermes_nemisis._data_directory) / "NEM_sci_packet_def.csv"
    )
    data = pkt.load(data_filename)
    return data


def l0_sci_data_to_l1(data: dict, original_filename: Path) -> Path:
    """
    Write level 0 nemisis science data to a level 1 cdf file.

    Parameters
    ----------
    data: dict
        A dictionary of arrays which includes the ccsds header fields
    original_filename: `pathlib.Path`
        The Path to the originating file.

    Returns
    -------
    output_filename: `pathlib.Path`
        Fully specificied filename of cdf file

    Examples
    --------
    >>> from pathlib import Path
    >>> from hermes_core.util.util import parse_science_filename
    >>> import hermes_nemisis.calibration as calib
    >>> data_filename = Path("hermes_NEM_l0_2022339-000000_v0.bin")
    >>> metadata = parse_science_filename(data_filename)  # doctest: +SKIP
    >>> data_packets = calib.parse_l0_sci_packets(data_filename)  # doctest: +SKIP
    >>> cdf_filename = calib.l0_sci_data_to_l1(data_packets, data_filename)  # doctest: +SKIP
    """
    ####  
    #### KRB: note that this code from the template is not used.
    ####      the file name is created by HermesData.save().
    ####      What we want to do is figure out the version number and 
    ####      time stamp for the L1 data file, IF they depends on the
    ####      L0 version number.
    # file_metadata = parse_science_filename(original_filename.name)

    # cdf_filename = original_filename.parent / create_science_filename(
    #     file_metadata["instrument"],
    #     file_metadata["time"],
    #     "l1",
    #     f'1.0.{file_metadata["version"]}',
    # )
    ######

    from astropy.time import Time, TimeDelta
    import astropy.units as u
    from astropy.timeseries import TimeSeries
    import numpy as np

    pckt_times = []
    for pckt in range(len(data['MET_SECONDS'])):
        pckt_times.append(data['MET_SECONDS'][pckt]+data['SCTF_SECONDS'][pckt]+data['SCTF_SUBSECONDS'][pckt]/2**32)

    pckt_delta = np.median(np.diff(pckt_times))
    fg_times = []
    fg_bx = []
    fg_by = []
    fg_bz = []
    pni1_bx = []
    pni1_by = []
    pni1_bz = []
    pni2_bx = []
    pni2_by = []
    pni2_bz = []
    for pckt in range(len(pckt_times)):
        for i in range(40):
            sample_time = pckt_times[pckt] + pckt_delta * (i - data['ACQ_POS_NUM'][pckt]) / 40 - data['ISR_COUNT'][pckt]*14.933e-6
            fg_times.append(sample_time)
            fg_bx.append(data['MAG_DATA'][pckt][i*9 ])
            fg_by.append(data['MAG_DATA'][pckt][i*9 +1])
            fg_bz.append(data['MAG_DATA'][pckt][i*9 +2])

            pni1_bx.append(data['MAG_DATA'][pckt][i*9 +3])
            pni1_by.append(data['MAG_DATA'][pckt][i*9 +4])
            pni1_bz.append(data['MAG_DATA'][pckt][i*9 +5])

            pni2_bx.append(data['MAG_DATA'][pckt][i*9 +6])
            pni2_by.append(data['MAG_DATA'][pckt][i*9 +7])
            pni2_bz.append(data['MAG_DATA'][pckt][i*9 +8])
    
    # convert fluxgate to pseudo-nT
    fg_bx = np.array(fg_bx) / 8388608.0 * 65000.0     
    fg_by = np.array(fg_by) / 8388608.0 * 65000.0
    fg_bz = np.array(fg_bz) / 8388608.0 * 65000.0

    # convert PNI to pseudo-nT
    pni1_bx = np.array(pni1_bx) * 13.66 / 7.0            
    pni1_by = np.array(pni1_by) * 13.66 / 7.0
    pni1_bz = np.array(pni1_bz) * 13.66 / 7.0
    pni2_bx = np.array(pni2_bx) * 13.66 / 7.0            
    pni2_by = np.array(pni2_by) * 13.66 / 7.0
    pni2_bz = np.array(pni2_bz) * 13.66 / 7.0

    # temperature conversions to engineering units
    t_eb = np.array(data['EBOX_TEMP']) / 128 
    t_fg = np.array(data['FG_TEMP']) * (3.3 / 1024) * 100 - 273.15
    t_pni1 = np.array(data['PNI1_TEMP']) # PNI MI-S1 sensor temperature conversion **GOES HERE**
    t_pni2 = np.array(data['PNI2_TEMP']) # PNI MI-S2 sensor temperature conversion **GOES HERE**

    # convert time tags to AstroPy times
    tt2k_epoch = Time('2000-01-01 11:58:55.816', scale='utc')
    pckt_time = tt2k_epoch + TimeDelta(pckt_times * u.s)
    fg_time = tt2k_epoch + TimeDelta(fg_times * u.s)

    # store 10 Sample/Sec (magnetometer cadence) data in AstroPy TimeSeries structures
    magts = TimeSeries(
      time=fg_time,
      data={
         'hermes_nem_fg_b1': u.Quantity(
            value = fg_bx,
            unit='nT',  # actually, 'pseudo-nT'  perhaps use u.def_unit('pseudo-nT') and u.add_enabled_units
            dtype=np.float32
         ),
         'hermes_nem_fg_b2': u.Quantity(
            value = fg_by,
            unit='nT',
            dtype=np.float32
         ),
         'hermes_nem_fg_b3': u.Quantity(
            value = fg_bz,
            unit='nT',
            dtype=np.float32
         ),

         # the PNIs should have their own timeseries, with its own time tags
         # otherwise, how do we account for the offset between the integration times?
         'hermes_nem_pni1_b1': u.Quantity(
            value = pni1_bx,
            unit='nT',  # actually, 'pseudo-nT'  perhaps use u.def_unit('pseudo-nT') and u.add_enabled_units
            dtype=np.float32
         ),
         'hermes_nem_pni1_b2': u.Quantity(
            value = pni1_by,
            unit='nT',
            dtype=np.float32
         ),
         'hermes_nem_pni1_b3': u.Quantity(
            value = pni1_bz,
            unit='nT',
            dtype=np.float32
         ),

         'hermes_nem_pni2_b1': u.Quantity(
            value = pni2_bx,
            unit='nT',  # actually, 'pseudo-nT'  perhaps use u.def_unit('pseudo-nT') and u.add_enabled_units
            dtype=np.float32
         ),
         'hermes_nem_pni2_b2': u.Quantity(
            value = pni2_by,
            unit='nT',
            dtype=np.float32
         ),
         'hermes_nem_pni2_b3': u.Quantity(
            value = pni2_bz,
            unit='nT',
            dtype=np.float32
         )

      }
   )

    # store 1 Sample/4 Sec (packet cadence) data in AstroPy TimeSeries structures
    pktts = TimeSeries(
        time=pckt_time,
        data={
            'hermes_nem_ebox_temp': u.Quantity(
            value = t_eb,
            unit='deg C',
            dtype=np.float32
            ),
            'hermes_nem_fg_temp': u.Quantity(
            value = t_fg,
            unit='deg C',
            dtype=np.float32
            ),
            'hermes_nem_pni1_temp': u.Quantity(
            value = t_pni1,
            unit='deg C',
            dtype=np.float32
            ),
            'hermes_nem_pni2_temp': u.Quantity(
            value = t_pni2,
            unit='deg C',
            dtype=np.float32
            )
        }
    )

    from hermes_core.timedata import HermesData
    from collections import OrderedDict

    ts = {
        'Epoch': magts,
        'Epoch_pkt': pktts
    }

    # now flesh out the HermesData structure
    input_attrs = HermesData.global_attribute_template("nemisis", "l1", "1.0.0")

    nemisis_data = HermesData(
        timeseries=ts,
        meta = input_attrs
    )

    nemisis_data.timeseries["Epoch"]["hermes_nem_fg_b1"].meta.update(
        OrderedDict({"CATDESC":"FG axis 1 OUTPUT IN PSEUDO-nT"})
    )
    nemisis_data.timeseries["Epoch"]["hermes_nem_fg_b2"].meta.update(
        OrderedDict({"CATDESC":"FG axis 2 OUTPUT IN PSEUDO-nT"})
    )
    nemisis_data.timeseries["Epoch"]["hermes_nem_fg_b3"].meta.update(
        OrderedDict({"CATDESC":"FG axis 3 OUTPUT IN PSEUDO-nT"})
    )

    nemisis_data.timeseries["Epoch"]["hermes_nem_pni1_b1"].meta.update(
        OrderedDict({"CATDESC":"PNI1 axis 1 OUTPUT IN PSEUDO-nT"})
    )
    nemisis_data.timeseries["Epoch"]["hermes_nem_pni1_b2"].meta.update(
        OrderedDict({"CATDESC":"PNI1 axis 2 OUTPUT IN PSEUDO-nT"})
    )
    nemisis_data.timeseries["Epoch"]["hermes_nem_pni1_b3"].meta.update(
        OrderedDict({"CATDESC":"PNI1 axis 3 OUTPUT IN PSEUDO-nT"})
    )

    nemisis_data.timeseries["Epoch"]["hermes_nem_pni2_b1"].meta.update(
        OrderedDict({"CATDESC":"PNI2 axis 1 OUTPUT IN PSEUDO-nT"})
    )
    nemisis_data.timeseries["Epoch"]["hermes_nem_pni2_b2"].meta.update(
        OrderedDict({"CATDESC":"PNI2 axis 2 OUTPUT IN PSEUDO-nT"})
    )
    nemisis_data.timeseries["Epoch"]["hermes_nem_pni2_b3"].meta.update(
        OrderedDict({"CATDESC":"PNI2 axis 3 OUTPUT IN PSEUDO-nT"})
    )


    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_ebox_temp"].meta.update(
        OrderedDict({"CATDESC":"EBOX temperature in deg C"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_fg_temp"].meta.update(
        OrderedDict({"CATDESC":"FG sensor temperature in deg C"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_pni1_temp"].meta.update(
        OrderedDict({"CATDESC":"PNI1 sensor temperature in deg C"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_pni2_temp"].meta.update(
        OrderedDict({"CATDESC":"PNI2 sensor temperature in deg C"})
    )


    cdf_file_path = nemisis_data.save(output_path=hermes_nemisis._data_directory, overwrite=True)
    return cdf_file_path
    # 
#    cdf_file_path = nemisis_data.save(output_path=cdf_filename, overwrite=True)


#    return cdf_filename


def get_calibration_file(data_filename: Path, time=None) -> Path:
    """
    Given a time, return the appropriate calibration file.

    Parameters
    ----------
    data_filename: `pathlib.Path`
        Fully specificied filename of the non-calibrated file (data level < 2)
    time: ~astropy.time.Time

    Returns
    -------
    calib_filename: `pathlib.Path`
        Fully specificied filename for the appropriate calibration file.

    Examples
    --------
    """
    return None


def read_calibration_file(calib_filename: Path):
    """
    Given a calibration, return the calibration structure.

    Parameters
    ----------
    calib_filename: `pathlib.Path`
        Fully specificied filename of the non-calibrated file (data level < 2)

    Returns
    -------
    output_filename: `pathlib.Path`
        Fully specificied filename of the appropriate calibration file.

    Examples
    --------
    """
    return None
