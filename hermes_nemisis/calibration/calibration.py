"""
A module for all things calibration.
"""

import random
from pathlib import Path

import ccsdspy
import numpy as np

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

    if (data_filename.stat().st_size == 0):
        raise ValueError(f"Zero length file: {data_filename}")
    data = pkt.load(data_filename)

    # calculate checksum of each packet
    n_pkts = len(data['CHECKSUM'])
    packet_bytes = ccsdspy.utils.split_packet_bytes(
        data_filename,
        include_primary_header=False
    )
    calculated_checksum_out = np.zeros(n_pkts, dtype=np.uint16)
    for i in range(n_pkts):
        CK_A = 0 
        CK_B = 0
        # each NEMISIS packet is 1117 bytest long.  
        # The last two bytes are the checksum bytes.
        # The checksum is calculatated on the first 1115 bytes.
        for octet_num in range(1115):
            CK_A = (CK_A + packet_bytes[i][octet_num+6]) % 255
            CK_B = (CK_B + CK_A) % 255
        calculated_checksum_out[i] = (CK_B << 8) | CK_A
    data["calculated_checksum"] = calculated_checksum_out

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
    log.info(f"Begin l0_sci_data_to_l1")

    # retrieve the version of the Level 0 input file, which 
    # consists of a single, two-digit version number.
    file_metadata = parse_science_filename(original_filename.name)
 
    # the Z-version of the Level 1 file  will reflect the version of the Level 0 file, 
    # but without any leading zero.
    cdf_version = f'1.0.{int(file_metadata["version"])}'

    from astropy.time import Time, TimeDelta
    import astropy.units as u
    from astropy.timeseries import TimeSeries

    n_pkts = len(data['MET_SECONDS'])
    acq_pos_times = np.zeros(n_pkts)
    acq_pos_num = np.zeros(n_pkts)
    max_pkt = 0

    # verify checksums
    checksum_valid = data["calculated_checksum"][range(n_pkts)] == data["CHECKSUM"][range(n_pkts)]

    ## Test bad packet(s) at the end of the file #### TEST
    ## poke in some data values to simulate bad packet(s)
    # data["START_FLAG"][n_pkts-1] = 1    #### TEST
    # data["XSUM"][n_pkts-2] = 0          #### TEST
    for pkt in range(n_pkts):
        # if RELAY_STATE == 0  OR  START_FLAG == 1  OR  XSUM == 0 then 
        # we must not use the the time-at-tone, ISR count or ACQ_POS_NUM.  
        if (data['RELAY_STATE'][pkt] == 1 and data['START_FLAG'][pkt] == 0 and data['XSUM'][pkt] == 1):
            # none of the above-mentioned conditions are true, so we trust the timing information
            # contained within this packet.
            # Record the precise time of the ACQ_POS_NUM sample of this packet
            # (we know the time at tone pulse occurs ISR_COUNT interrupts after the ACQ_POS_NUM'th sample is begun)
            acq_pos_times[pkt] = data['MET_SECONDS'][pkt] + data['SCTF_SECONDS'][pkt] + (data['SCTF_SUBSECONDS'][pkt]/2**32) - \
                (data['ISR_COUNT'][pkt]*14.933e-6)
            acq_pos_num[pkt] = data['ACQ_POS_NUM'][pkt]
            max_pkt = pkt
        else:
            log.info(f"ignoring packet time for PKT_ID {data['PKT_ID'][pkt]}, RELAY_STATE: {data['RELAY_STATE'][pkt]}, START_FLAG: {data['START_FLAG'][pkt]}, XSUM: {data['XSUM'][pkt]}")

    samp_delta = np.median(np.diff(acq_pos_times))/40

    # now find start-up packets and fill them in with time tags calculated from the following packet, 
    # working backwards from the last valid packet
    for pkt in range(max_pkt-1,-1,-1):
        if (acq_pos_times[pkt] == 0):
            acq_pos_times[pkt] = acq_pos_times[pkt+1] - samp_delta*40
            acq_pos_num[pkt] = acq_pos_num[pkt+1]

    # if the last one is bad, then we have to toss it:
    if (max_pkt + 1 < n_pkts):
        log.info(f"Truncating {n_pkts - max_pkt - 1} packets from end of file: PKT_IDs {data['PKT_ID'][range(max_pkt+1,n_pkts)]}")
        n_pkts = max_pkt + 1
        acq_pos_times = acq_pos_times[range(n_pkts)]
        acq_pos_num = acq_pos_num[range(n_pkts)]

    fg_times = np.zeros(n_pkts*40)
    fg_bx =    np.zeros(n_pkts*40)
    fg_by =    np.zeros(n_pkts*40)
    fg_bz =    np.zeros(n_pkts*40)
    pni1_bx =  np.zeros(n_pkts*40)
    pni1_by =  np.zeros(n_pkts*40)
    pni1_bz =  np.zeros(n_pkts*40)
    pni2_bx =  np.zeros(n_pkts*40)
    pni2_by =  np.zeros(n_pkts*40)
    pni2_bz =  np.zeros(n_pkts*40)
    sampflags =np.zeros(n_pkts*40, dtype = np.int8) 
    for pkt in range(n_pkts):
        # with bit 0 of the sampflags byte,
        # flag the first sample of the packet if the "power-on start" flag is present
        sampflags[pkt*40] |= data['START_FLAG'][pkt] << 0
        for i in range(40):
            sample_time = acq_pos_times[pkt] + samp_delta * (i - acq_pos_num[pkt])
            fg_times[pkt*40+i] = sample_time
            fg_bx[pkt*40+i] = data['MAG_DATA'][pkt][i*9 ]
            fg_by[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +1]
            fg_bz[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +2]

            pni1_bx[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +3]
            pni1_by[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +4]
            pni1_bz[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +5]

            pni2_bx[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +6]
            pni2_by[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +7]
            pni2_bz[pkt*40+i] = data['MAG_DATA'][pkt][i*9 +8]

            # with bit 1 of the sampflags byte,
            # flag all samples in the packet if the checksum is bad
            sampflags[pkt*40+i] |= (checksum_valid[pkt] == 0) << 1


    # end time of each packet (in theory)
    pkt_times = acq_pos_times + (40 - acq_pos_num)*samp_delta


    # convert fluxgate to pseudo-nT
    fg_bx = fg_bx / 8388608.0 * 65000.0     
    fg_by = fg_by / 8388608.0 * 65000.0
    fg_bz = fg_bz / 8388608.0 * 65000.0

    # convert PNI to pseudo-nT
    pni1_bx = pni1_bx * 13.66 / 7.0            
    pni1_by = pni1_by * 13.66 / 7.0
    pni1_bz = pni1_bz * 13.66 / 7.0
    pni2_bx = pni2_bx * 13.66 / 7.0            
    pni2_by = pni2_by * 13.66 / 7.0
    pni2_bz = pni2_bz * 13.66 / 7.0

    # temperature conversions to engineering units
    t_eb = np.array(data['EBOX_TEMP'][range(n_pkts)]) / 128 
    t_fg = np.array(data['FG_TEMP'][range(n_pkts)]) * (3.3 / 1024) * 100 - 273.15
    t_pni1 = np.array(data['PNI1_TEMP'][range(n_pkts)]) # PNI MI-S1 sensor temperature conversion **GOES HERE**
    t_pni2 = np.array(data['PNI2_TEMP'][range(n_pkts)]) # PNI MI-S2 sensor temperature conversion **GOES HERE**

    # convert time tags to AstroPy times
    tt2k_epoch = Time('2000-01-01 11:58:55.816', scale='utc')
    pkt_time = tt2k_epoch + TimeDelta(pkt_times * u.s)
    fg_time = tt2k_epoch + TimeDelta(fg_times * u.s)

    # RT clock times
    rtc_year = bcd_to_int(np.array(data["RTC_YEAR"][range(n_pkts)]))  #  convert from Binary Coded Decimal
    rtc_month = bcd_to_int(np.array(data["RTC_MONTH"][range(n_pkts)]))
    rtc_day = bcd_to_int(np.array(data["RTC_DAY"][range(n_pkts)]))
    rtc_hours = bcd_to_int(np.array(data["RTC_HOURS"][range(n_pkts)]))
    rtc_minutes = bcd_to_int(np.array(data["RTC_MINUTES"][range(n_pkts)]))
    rtc_seconds = bcd_to_int(np.array(data["RTC_SECONDS"][range(n_pkts)]))
    rt_clock = []
    for pkt in range(n_pkts):
        rt_clock_time = Time(
            "20%02d-%02d-%02d %02d:%02d:%02d.000" % (
                rtc_year[pkt], rtc_month[pkt], rtc_day[pkt], 
                rtc_hours[pkt], rtc_minutes[pkt], rtc_seconds[pkt]
            ), 
            scale='utc'
        )
        rt_clock.append(rt_clock_time.gps)

    log.info(f"Creating 10 S/s TimeSeries structure")

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
         ),

         'hermes_nem_sampflags': u.Quantity(
            value = sampflags,
            unit='',
            dtype=np.uint8
         )

      }
   )

    log.info(f"Creating 1/4 S/s TimeSeries structure")

    # store 1 Sample/4 Sec (packet cadence) data in AstroPy TimeSeries structures
    pktts = TimeSeries(
        time=pkt_time,
        data={
            'hermes_nem_pkt_id': u.Quantity(
                value = data['PKT_ID'][range(n_pkts)],
                unit='',
                dtype=np.uint32
            ),
            'hermes_nem_acq_pos_num': u.Quantity(
                value = data['ACQ_POS_NUM'][range(n_pkts)],
                unit='',
                dtype=np.uint8
            ),
            'hermes_nem_xsum': u.Quantity(
                value = data['XSUM'][range(n_pkts)],
                unit='',
                dtype=np.uint8
            ),
            'hermes_nem_start_flag': u.Quantity(
                value = data['START_FLAG'][range(n_pkts)],
                unit='',
                dtype=np.uint8
            ),
            'hermes_nem_ops_mode': u.Quantity(
                value = data['OPS_MODE'][range(n_pkts)],
                unit='',
                dtype=np.uint8
            ),
            'hermes_nem_relay_state': u.Quantity(
                value = data['RELAY_STATE'][range(n_pkts)],
                unit='',
                dtype=np.uint8
            ),
            'hermes_nem_rtc_seconds': u.Quantity(
                value = rtc_seconds,
                unit='s',
                dtype=np.uint16
            ),
            'hermes_nem_rtc_minutes': u.Quantity(
                value = rtc_minutes,
                unit='min',
                dtype=np.uint16
            ),
            'hermes_nem_rtc_hours': u.Quantity(
                value = rtc_hours,
                unit='hr',
                dtype=np.uint16
            ),
            'hermes_nem_rtc_day': u.Quantity(
                value = rtc_day,
                unit='day',
                dtype=np.uint16
            ),
            'hermes_nem_rtc_month': u.Quantity(
                value = rtc_month,
                unit='',
                dtype=np.uint16
            ),
            'hermes_nem_rtc_year': u.Quantity(
                value = rtc_year,
                unit='',
                dtype=np.uint16
            ),
            'hermes_nem_rtc_time': u.Quantity(
                value = rt_clock,
                unit='s',
                dtype=np.float64
            ),
            'hermes_nem_met_seconds': u.Quantity(
                value = data['MET_SECONDS'][range(n_pkts)],
                unit='s',
                dtype=np.uint32
            ),
            'hermes_nem_sctf_seconds': u.Quantity(
                value = data['SCTF_SECONDS'][range(n_pkts)],
                unit='s',
                dtype=np.uint32
            ),
            'hermes_nem_sctf_subseconds': u.Quantity(
                value = data['SCTF_SUBSECONDS'][range(n_pkts)]/2**32,
                unit='s',
                dtype=np.float64
            ),
            'hermes_nem_isr_count': u.Quantity(
                value = data["ISR_COUNT"][range(n_pkts)],
                unit='',
                dtype=np.uint16
            ),
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
            ),
            'hermes_nem_checksum_valid': u.Quantity(
                value = checksum_valid,
                unit='',
                dtype=np.uint8
            ),
            'hermes_nem_checksum': u.Quantity(
                value = data["CHECKSUM"][range(n_pkts)],
                unit='',
                dtype=np.uint16
            ),
            'hermes_nem_calculated_checksum': u.Quantity(
                value = data["calculated_checksum"][range(n_pkts)],
                unit='',
                dtype=np.uint16
            )
        }
    )

    log.info(f"Creating HermesData structure")

    from hermes_core.timedata import HermesData
    from collections import OrderedDict

    ts = {
        'Epoch': magts,
        'Epoch_pkt': pktts
    }

    # now flesh out the HermesData structure
    input_attrs = HermesData.global_attribute_template("nemisis", "l1", cdf_version)

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

    nemisis_data.timeseries["Epoch"]["hermes_nem_sampflags"].meta.update(
        OrderedDict({"CATDESC":"FLAGS: bit 0: power-on, bit 1: bad cksum"})
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

    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_pkt_id"].meta.update(
        OrderedDict({"CATDESC":"packet ID"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_acq_pos_num"].meta.update(
        OrderedDict({"CATDESC":"Acquisition position number"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_xsum"].meta.update(
        OrderedDict({"CATDESC":"MOSI Packet Checksum Comparison Agreed"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_start_flag"].meta.update(
        OrderedDict({"CATDESC":"Power-On Start Flag"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_ops_mode"].meta.update(
        OrderedDict({"CATDESC":"Last Packet before shutdown"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_relay_state"].meta.update(
        OrderedDict({"CATDESC":"Transorb Relay Engaged"})
    )

    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_rtc_seconds"].meta.update(
        OrderedDict({"CATDESC":"RTC seconds"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_rtc_minutes"].meta.update(
        OrderedDict({"CATDESC":"RTC minutes"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_rtc_hours"].meta.update(
        OrderedDict({"CATDESC":"RTC hours"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_rtc_day"].meta.update(
        OrderedDict({"CATDESC":"RTC day"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_rtc_month"].meta.update(
        OrderedDict({"CATDESC":"RTC month"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_rtc_year"].meta.update(
        OrderedDict({"CATDESC":"RTC year"})
    )

    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_met_seconds"].meta.update(
        OrderedDict({"CATDESC":"MET seconds"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_sctf_seconds"].meta.update(
        OrderedDict({"CATDESC":"SCTF seconds"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_sctf_subseconds"].meta.update(
        OrderedDict({"CATDESC":"SCTF subseconds"})
    )

    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_rtc_time"].meta.update(
        OrderedDict({"CATDESC":"RTC time as a GPS time"})
    )

    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_isr_count"].meta.update(
        OrderedDict({"CATDESC":"Interrupt service request count"})
    )

    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_checksum_valid"].meta.update(
        OrderedDict({"CATDESC":"telemetry packet checksum valid"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_checksum"].meta.update(
        OrderedDict({"CATDESC":"telemetry packet checksum"})
    )
    nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_calculated_checksum"].meta.update(
        OrderedDict({"CATDESC":"telemetry packet calculated checksum"})
    )
    log.info(f"Finished Creating HermesData structure")

    cdf_file_path = nemisis_data.save(output_path=hermes_nemisis._data_directory, overwrite=True)

    log.info(f"Wrote l1 cdf file:{cdf_file_path.name}.")

    return cdf_file_path

def bcd_to_int(bcd):
    """
    Given a byte in Binary Coded Decimal, return the decoded value
    will work for a scalar input or a numpy array of input values
    """

    return (bcd >> 4) * 10 + (bcd & 0xF)

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
