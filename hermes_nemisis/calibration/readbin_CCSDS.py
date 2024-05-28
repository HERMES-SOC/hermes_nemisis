def readbin_CCSDS(binfile):
   import numpy as np
   import ccsdspy

   #************************************************************************************
   #
   #                                    R E A D B I N
   #
   #  Program:      READBIN - for plain CCSDS file
   #
   #  Programmer:   David G. Simpson & Kenneth R Bromund
   #                NASA Goddard Space Flight Center
   #                Greenbelt, Maryland  20771
   #
   #  Date:         June 2, 2023  
   #                Modified by KRB May 15, 2024 to use HermesData to write to CDF
   #
   #  Language:     Python 3
   #
   #  Version:      2.0
   #
   #  Description:  Read binary HERMES data files in plain CCSDS format.
   #
   #                The program loops over the bytes in one packet, converts them to
   #                ints, and stores them in a list.  Once an entire packet is
   #                collected, the list is decoded.  We then repeat with the next
   #                packet.
   #
   #                It is not easy to determine data sizes in Python.  For that
   #                reason, we need a trick when dealing with negative 2's complement
   #                integers:  if bit 23 = 1, then replace the integer A with
   #                   -(2^24 - A)
   #
   #************************************************************************************

   # Set of initial lists which will be appended with data from each good packet
   fg_bx = []  
   fg_by = []
   fg_bz = []
   pni1_bx = []
   pni1_by = []
   pni1_bz = []
   pni2_bx = []
   pni2_by = []
   pni2_bz = []
   attone_met_out = []
   attone_stcf_sec_out = []
   attone_stcf_subsec_out = []
   acqpos_out = []
   isrcount_out = []
   sync_out = []
   ebox_temp_out = []
   fg_temp_out = []
   pni1_temp_out = []
   pni2_temp_out = []
   cksum_out = []
   calculated_cksum_out = []

   # Set up dictionaries to store data in a per-packet format; can be useful when seeing if specific samples from each packet are corrupted.
   fgx_per_packet = dict.fromkeys(np.arange(40))
   fgy_per_packet = dict.fromkeys(np.arange(40))
   fgz_per_packet = dict.fromkeys(np.arange(40))
   for i, _ in fgx_per_packet.items():
      fgx_per_packet[i] = []
      fgy_per_packet[i] = []
      fgz_per_packet[i] = []
   

   #
   #  After skipping past initial "mystery bytes", we read the packet contents
   #  byte-by-byte into a list called blist[].  Once all bytes for a packet
   #  have been read in, we decode all the bytes in blist[], then move on
   #  to the next packet and repeat.
   #
   n_bad_sync = 0
   n_bad_cksum = 0
   total_pnum = 0
   # pkt = ccsdspy.FixedLength.from_file(
   #    os.path.join(hermes_nemisis._data_directory, "hermes_nemisis/data/MAG_sci_packet_def.csv")
   # )
   # data = pkt.load(data_filename)

   with open(binfile, "rb") as f:                 # loop on bytes in binary file
      bnum = -1                                   # initialize byte number
      lctr = 0                                    # initialie loop counter
      pnum = 0                                    # init packet number
      blist = []                                  # init list of bytes
      while (byte := f.read(1)):                  # get byte from file
         bnum += 1                                # increment byte counter
         if (bnum < 12):                          # skip past CCSDS HEADER and time tag
            continue
         ibyte = int.from_bytes(byte, 'big')      # convert byte to int
         blist.append(ibyte)                      # append int to end of blist[]
         lctr += 1                                # increment loop counter
         if (lctr == (1117+12)):                  # here we've read  1117 byte NEMISIS packet + 6 byte ccsds header + 6 byte time for the next packet
                                                   #  blist[], so start decoding
            pnum += 1                             # increment packet number
            total_pnum += 1                       # increment total number of packets processed

            sync = (blist[0] << 8) | blist[1]     # construct sync word
            sync_out.append(sync)

            cksum = (blist[1115] << 8) | blist[1116]   # construct checksum
            cksum_out.append(cksum)

            # Calculate checksum
            CK_A = 0 
            CK_B = 0
            for octet_num in range(1115):
               CK_A = (CK_A + blist[octet_num]) % 255
               CK_B = (CK_B + CK_A) % 255
            calculated_cksum = (CK_B << 8) | CK_A
            calculated_cksum_out.append(calculated_cksum)
            
            # check for correct sync; if sync is bad, no data will be appended to the lists except for the calculated checksum, provided checksum, and sync.
            if (sync != 59379):
               n_bad_sync += 1
               print('Bad Sync at Packet Num %d' % total_pnum)
            # check for correct checksum; if checksum is bad, no data will be appended to the lists except for the calculated checksum, provided checksum, and sync.
            elif (calculated_cksum != cksum):
               n_bad_cksum += 1
               print('Bad Checksum at Packet Num %d' % total_pnum)
            # if checksum and sync are both good, start appending data from the packet.
            elif (calculated_cksum == cksum) and (sync == 59379):
               pid = blist[2]                        # construct PID
               pid = (pid << 8) | blist[3]
               pid = (pid << 8) | blist[4]
               pid = (pid << 8) | blist[5]

               acqpos = blist[6]                     # get ACQPOS
               acqpos_out.append(acqpos)

               opstat = blist[7]                     # get OPSTAT

               rtc_csec = blist[8]                   # get RTC clock bytes
               rtc_sec = blist[9]
               rtc_min = blist[10]
               rtc_hr = blist[11]
               rtc_day = blist[12]
               rtc_mon = blist[13]
               rtc_year = blist[14]

               attone_met = blist[15]                      # construct AtTone MET clock
               attone_met = (attone_met << 8) | blist[16]
               attone_met = (attone_met << 8) | blist[17]
               attone_met = (attone_met << 8) | blist[18]
               attone_met_out.append(attone_met)

               attone_stcf_sec = blist[19]                 # construct AtTone STCF seconds
               attone_stcf_sec = (attone_stcf_sec << 8) | blist[20]
               attone_stcf_sec = (attone_stcf_sec << 8) | blist[21]
               attone_stcf_sec = (attone_stcf_sec << 8) | blist[22]

               attone_stcf_ss = blist[23]                          # construct AtTone STCF subsec
               attone_stcf_ss = (attone_stcf_ss << 8) | blist[24]
               attone_stcf_ss = (attone_stcf_ss << 16)
               attone_stcf_sec_out.append(attone_stcf_sec)
               attone_stcf_subsec_out.append(attone_stcf_ss)

               isrcount = blist[25]                       # construct isr count
               isrcount = (isrcount << 8) | blist[26]
               isrcount_out.append(isrcount)

               for i in range(40):                        # loop on magnetometer data
                  bfg_x = blist[27*i+27]                  # construct fluxgate Bx
                  bfg_x = (bfg_x << 8) | blist[27*i+28]
                  bfg_x = (bfg_x << 8) | blist[27*i+29]
                  if ((bfg_x & 0x800000) != 0):           # if negative..
                     bfg_x = -(0x1000000 - bfg_x)         # 2's compl
                  bfg_y = blist[27*i+30]                  # construct fluxgate By
                  bfg_y = (bfg_y << 8) | blist[27*i+31]
                  bfg_y = (bfg_y << 8) | blist[27*i+32]
                  if ((bfg_y & 0x800000) != 0):           # if negative..
                     bfg_y = -(0x1000000 - bfg_y)         # 2's compl
                  bfg_z = blist[27*i+33]                  # construct fluxgate Bz
                  bfg_z = (bfg_z << 8) | blist[27*i+34]
                  bfg_z = (bfg_z << 8) | blist[27*i+35]
                  if ((bfg_z & 0x800000) != 0):           # if negative..
                     bfg_z = -(0x1000000 - bfg_z)         # 2's compl

                  bfg_x = bfg_x / 8388608.0 * 65000.0     # apply fluxgate calibration curves
                  bfg_y = bfg_y / 8388608.0 * 65000.0
                  bfg_z = bfg_z / 8388608.0 * 65000.0
                  fg_bx.append(bfg_x)
                  fg_by.append(bfg_y)
                  fg_bz.append(bfg_z)
                  fgx_per_packet[i].append(bfg_x)
                  fgy_per_packet[i].append(bfg_y)
                  fgz_per_packet[i].append(bfg_z)

                  bmis1_x = blist[27*i+36]                   # construct PNI MI-S1 Bx
                  bmis1_x = (bmis1_x << 8) | blist[27*i+37]
                  bmis1_x = (bmis1_x << 8) | blist[27*i+38]
                  if ((bmis1_x & 0x800000) != 0):            # if negative..
                     bmis1_x = -(0x1000000 - bmis1_x)        # 2's compl
                  bmis1_y = blist[27*i+39]                   # construct PNI MI-S1 By
                  bmis1_y = (bmis1_y << 8) | blist[27*i+40]
                  bmis1_y = (bmis1_y << 8) | blist[27*i+41]
                  if ((bmis1_y & 0x800000) != 0):            # if negative..
                     bmis1_y = -(0x1000000 - bmis1_y)        # 2's compl
                  bmis1_z = blist[27*i+42]                   # construct PNI MI-S1 Bz
                  bmis1_z = (bmis1_z << 8) | blist[27*i+43]
                  bmis1_z = (bmis1_z << 8) | blist[27*i+44]
                  if ((bmis1_z & 0x800000) != 0):            # if negative..
                     bmis1_z = -(0x1000000 - bmis1_z)        # 2's compl

                  bmis1_x = bmis1_x * 13.66 / 7.0            # ** apply TEMPORARY PNI calibration curves
                  bmis1_y = bmis1_y * 13.66 / 7.0
                  bmis1_z = bmis1_z * 13.66 / 7.0
                  pni1_bx.append(bmis1_x)
                  pni1_by.append(bmis1_y)
                  pni1_bz.append(bmis1_z)

                  bmis2_x = blist[27*i+45]                   # construct PNI MI-S2 Bx
                  bmis2_x = (bmis2_x << 8) | blist[27*i+46]
                  bmis2_x = (bmis2_x << 8) | blist[27*i+47]
                  if ((bmis2_x & 0x800000) != 0):            # if negative..
                     bmis2_x = -(0x1000000 - bmis2_x)        # 2's compl
                  bmis2_y = blist[27*i+48]                   # construct PNI MI-S2 By
                  bmis2_y = (bmis2_y << 8) | blist[27*i+49]
                  bmis2_y = (bmis2_y << 8) | blist[27*i+50]
                  if ((bmis2_y & 0x800000) != 0):            # if negative..
                     bmis2_y = -(0x1000000 - bmis2_y)        # 2's compl
                  bmis2_z = blist[27*i+51]                   # construct PNI MI-S2 Bz
                  bmis2_z = (bmis2_z << 8) | blist[27*i+52]
                  bmis2_z = (bmis2_z << 8) | blist[27*i+53]
                  if ((bmis2_z & 0x800000) != 0):            # if negative..
                     bmis2_z = -(0x1000000 - bmis2_z)        # 2's compl

                  bmis2_x = bmis2_x * 13.66 / 7.0            # ** apply TEMPORARY PNI calibration curves
                  bmis2_y = bmis2_y * 13.66 / 7.0
                  bmis2_z = bmis2_z * 13.66 / 7.0
                  pni2_bx.append(bmis2_x)
                  pni2_by.append(bmis2_y)
                  pni2_bz.append(bmis2_z)

               t_eb = blist[1107]                    # construct electronics box temperature
               t_eb = (t_eb << 8) | blist[1108]
               if ((t_eb & 0x8000) != 0):            # if negative..
                  t_eb = -(0x10000 - t_eb)           # 2's compl
               t_eb = t_eb / 128                     
               ebox_temp_out.append(t_eb)

               t_fg = blist[1109]                    # construct fluxgate sensor temperature
               t_fg = (t_fg << 8) | blist[1110]
               if ((t_fg & 0x8000) != 0):            # if negative..
                  t_fg = -(0x10000 - t_fg)           # 2's compl
               t_fg = t_fg * (3.3 / 1024) * 100 - 273.15                    
               fg_temp_out.append(t_fg)

               t_mis1 = blist[1111]                  # construct PNI MI-S1 sensor temperature
               t_mis1 = (t_mis1 << 8) | blist[1112]
               if ((t_mis1 & 0x8000) != 0):          # if negative..
                  t_mis1 = -(0x10000 - t_mis1)       # 2's compl
               pni1_temp_out.append(t_mis1)

               t_mis2 = blist[1113]                  # construct PNI MI-S2 sensor temperature
               t_mis2 = (t_mis2 << 8) | blist[1114]
               if ((t_mis2 & 0x8000) != 0):          # if negative..
                  t_mis2 = -(0x10000 - t_mis2)       # 2's compl
               pni2_temp_out.append(t_mis2)

            lctr = 0                              # re-initialize loop counter
            blist = []                            # re-initialize list
   
   from astropy.time import Time, TimeDelta
   import astropy.units as u
   from astropy.timeseries import TimeSeries

   pckt_times = []
   tt2k_epoch = Time('2000-01-01 11:58:55.816', scale='utc')
   for pckt in range(len(attone_met_out)):
      pckt_times.append(attone_met_out[pckt]+attone_stcf_sec_out[pckt]+attone_stcf_subsec_out[pckt]/2**32)

   pckt_delta = np.median(np.diff(pckt_times))
   fg_times = []
   for pckt in range(len(pckt_times)):
      for i in range(40):
         sample_time = pckt_times[pckt] + pckt_delta * (i - acqpos_out[pckt]) / 40 - isrcount_out[pckt]*14.933e-6
         fg_times.append(sample_time)

   pckt_time = tt2k_epoch + TimeDelta(pckt_times * u.s)
   fg_time = tt2k_epoch + TimeDelta(fg_times * u.s)

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


   tempts = TimeSeries(
      time=pckt_time,
      data={
         'hermes_nem_ebox_temp': u.Quantity(
            value = ebox_temp_out,
            unit='deg C',
            dtype=np.float32
         ),
         'hermes_nem_fg_temp': u.Quantity(
            value = fg_temp_out,
            unit='deg C',
            dtype=np.float32
         ),
         'hermes_nem_pni1_temp': u.Quantity(
            value = pni1_temp_out,
            unit='deg C',
            dtype=np.float32
         ),
         'hermes_nem_pni2_temp': u.Quantity(
            value = pni2_temp_out,
            unit='deg C',
            dtype=np.float32
         )
      }
   )

   from hermes_core.timedata import HermesData
   from collections import OrderedDict

   ts = {
      'Epoch': magts,
      'Epoch_pkt': tempts
   }

   input_attrs = HermesData.global_attribute_template("nemisis", "l1", "1.0.0")
   # here is where I test Andrew's new code to allow more than one timeseries
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
      OrderedDict({"CATDESC":"PNI2 axis 3 OUTPUT IN PSEUDO-nT"})
   )
   nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_fg_temp"].meta.update(
      OrderedDict({"CATDESC":"PNI2 axis 3 OUTPUT IN PSEUDO-nT"})
   )
   nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_pni1_temp"].meta.update(
      OrderedDict({"CATDESC":"PNI2 axis 3 OUTPUT IN PSEUDO-nT"})
   )
   nemisis_data.timeseries["Epoch_pkt"]["hermes_nem_pni2_temp"].meta.update(
      OrderedDict({"CATDESC":"PNI2 axis 3 OUTPUT IN PSEUDO-nT"})
   )


   cdf_file_path = nemisis_data.save(output_path="./", overwrite=True)

   return cdf_file_path