.. _data:

****
Data
****

Overview
========

Data Description
----------------

+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+
| Level | Product                                           | Description                                                                               |
+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+
| 1     | Vector Magnetic Field from 3 sensors (M0, M1, M2) | CCSDS, each packet is 4-sec long at 10 Hz rate, 3-axis magnetic field components          |
|       | in spacecraft coordinates                         | for all three sensors in coordinate system native to the HERMES Payload                   |
+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+
| 1     | Sensor Temperatures                               | Temperatures at M0, M1, M2                                                                |
+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+
| 2     | Vector Magnetic Field from 3 sensors (M0, M1, M2) | 3d vector magnetic fields (Bx, By, Bz) in nT for each magnetometer (M0, M1, M2)           |
|       | in GSE coordinates                                | using final calibrations for offsets and gains                                            |
+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+
| 2     | Sensor Temperatures                               | Final calibrated temperatures                                                             |
+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+
| 3     | Magnetic field at Gateway                         | Background-subtracted and processed 3d vector magnetic field (Bx, By, Bz) in nT in a      |
|       |                                                   | common coordinate system (e.g. GSE). Derived by combining individual sensor data.         |
+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+
| QL    | Vector Magnetic Field from 3 sensors (M0, M1, M2) | Despiked values for 3d vector magnetic fields (Bx, By, Bz) in nT for each magnetometer    |
|       |  in GSE coordinates (Unvalidated)                 | (M0, M1, M2) in their local coordinate system, plus temperatures in Celsius for each      |
|       |                                                   | sensor. (and time-corrected and time-checked)                                             |
+-------+---------------------------------------------------+-------------------------------------------------------------------------------------------+




Getting Data
============



Reading Data
============



Calibrating Data
================
Data products below level 2 generally require calibration to be transformed into scientificically useable units.
This section describes how to calibrate data files from lower to higher levels.