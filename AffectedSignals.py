# -*- coding: utf-8 -*-
"""
Created on Sat May 11 30:14:43 2020

@author: lauraigs and vedantm
"""

#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""

The following script is meant to tell the user about the usability of different signals due to interference of satellite signals.
The user needs start and end times in MJD, and a signal.
An example input would be:
Start Time (MJD): 58961.931192
Stop Time (MJD): 58961.933831
Signal ID: 23253
This needs improvement, especially when reading in different signals. There is a problem reading planetlist.txt because of the different formatting. All you would need is the RA and DEC of the different sources placed in the database. I assume there is a better way to find this. 

"""

# grabbing stuff from the database
from setiS20.sql import get_data_with_clause
from setiS20.sql import get_data
from setiS20 import CandidateSignal
from setiS20.globals import (
    FFT_LENGTH,
    TIME_PER_ROW,
    OBSERVATIONS,
)  # may take 5 seconds or more to execute - be patient!
import os
import sys
import spiceypy as spice

# this gets rid of the warnings! be very careful!
import warnings

warnings.filterwarnings("ignore")

sys.path.append(".")
from SatPos import TLEImportReturn
from SatPos import DopplerFactorReturn
import SatPos
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.utils.iers import conf

conf.iers_auto_url = "http://storage.googleapis.com/epss179s20/finals2000A.all"
planetlist = np.loadtxt("planetlist.txt", skiprows=7, dtype=str, unpack=True)


# functions:


def radec(ra, dec):
    """
    This returns right ascension and declination in radians
    """
    s = SkyCoord(ra + dec, unit=(u.hourangle, u.deg), frame="icrs")
    ra = s.ra.radian
    dec = s.dec.radian
    return ra, dec


def abs_separation(ra1, dec1, ra2, dec2):
    """
    Given the right ascension and declination of two coordinates, this function uses astropy to return the absolute separation.
    """
    c1 = SkyCoord(ra1 * u.rad, dec1 * u.rad, frame="icrs")
    c2 = SkyCoord(ra2 * u.rad, dec2 * u.rad, frame="icrs")
    sep = c1.separation(c2)
    return sep


def altitude(time, ra, dec):
    """
    This returns altitude with respect to GBT given right ascension, declination, and time
    """
    s = SkyCoord(ra=time * u.rad, dec=dec * u.rad, frame="icrs")
    time = Time(time, format="mjd", scale="utc")
    loc = EarthLocation.of_site("Green Bank Telescope")
    attribute = AltAz(obstime=time, location=loc)
    altaz = s.transform_to(attribute)
    alt = altaz.alt.value
    return alt


##Prompt
starttime = str(input("Start Time (MJD): "))
endtime = str(input("Stop Time (MJD): "))
signal_id = input("Signal ID: ")
timewindow = [starttime, endtime]
print("Name              Possible Interference   Separation     ")


# initial parameters
utctimewindow = Time(timewindow, format="mjd", scale="utc")
utcrange = utctimewindow.to_value("iso")
step = 1000
obs_long = -79.83983858
obs_lat = 38.24590630
f0 = 15e9  # IMPROVE ME!! what is the center frequency of the GPS satellite? How do we find this information????
obs_rad = 6370.7403920


# Obtaining multiple TLE files inside a folder.
names, tle_w = SatPos.m_tle_import(
    "TLEs", utcrange, step, obs_rad, obs_long, obs_lat, f0
)
for name, tlereturn in zip(names, tle_w):
    sv = tlereturn.sat_sv[:3]
    time = tlereturn.mjd_time

    # The following needs improvement. It is just collecting the source's information.
    s = CandidateSignal(signal_id)  # identify the signal
    sloc = [
        planetlist[0] == (s.source + "_" + str(s.scan))
    ]  # identify the source and the source's location
    source_ra, source_dec = planetlist[1][sloc], planetlist[2][sloc]  # ra/dec
    source_ra, source_dec = radec(source_ra[0], source_dec[0])  # ra/dec in radians

    # This loop finds each state vector at a given time. From there, it tells you the approximate separation.
    flags = 0
    for i in range(len(time)):
        mod_sv = [sv[0][i], sv[1][i], sv[2][i]]  # state vector from the satellite
        ra_dec = spice.spiceypy.recrad(
            [sv[0][i], sv[1][i], sv[2][i]]
        )  # obtaining ra/dec from the state vector using spiceypy
        alt = altitude(time[i], ra_dec[1], ra_dec[2])
        if alt > 0:
            sep = abs_separation(source_ra, source_dec, ra_dec[1], ra_dec[2])
            print("{0:5}      Yes                     {1:3.5} deg".format(name, sep))
            flags = 1
            break

    if flags == 0:
        print("{0:5}      No                     ".format(name))