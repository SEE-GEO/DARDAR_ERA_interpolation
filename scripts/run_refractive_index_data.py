#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:00:30 2020

Script to download ERA5 data for generating the refractive index data
for Ericsson. Generates pickle files containing the ERA5 parameters,
The conversion to csv files with refractive index is made with the script
write_csv_ericsson.py

@author: inderpreet
"""
import os
import glob
import numpy as np
import random
from era2dardar.ERA5 import ERA5p, ERA5s
#from era2dardar.dardar2atmdata import dardar2atmdata
from era2dardar.DARDAR import DARDARProduct
from era2dardar.locations import locations
from era2dardar.RADARLIDAR import DARDAR, CLOUDSAT
from era2dardar.utils.alt2pressure import alt2pres
#import typhon.arts.xml as xml
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import shutil
import pickle
from era2dardar.utils.match_dardar_cloudsat import match_dardar_cloudsat



def run_all_cases(p_grid, onsala, t0, outpath):


        outdir = t0.strftime("%Y") + "_" + t0.strftime("%3j") + "_" + t0.strftime("%2H")

        print (outdir)
        if os.path.isfile(os.path.join(outpath, outdir + '.zip')):
            print ('file %s already exists, doing next file'%outdir)
            t0 +=  timedelta(days = 1)
            return


        print ("t_0, t_1", onsala.t_0, onsala.t_1)


        lon          = onsala.longitude
        lat          = onsala.latitude

        lat1         = np.around(lat.min() - 1)
        lat2         = np.around(lat.max() + 1)

        lon1         = np.around(lon.min() - 1)
        lon2         = np.around(lon.max() + 1)

        if lon1 < -180:
            lon1 = -180.0
        if lon2 > 180:
            lon2 = 180.0

        # when encountering prime meridian, download global data
        # keeps interpolation simple
        if lon1 == -180 or lon2 == 180:
            lon1 = -180.
            lon2 =  180.

        # domain for which ERA5 data is downloaded
        domain  = [lat1, lat2, lon1, lon2]

        #eras = ERA5s(onsala.t_0, onsala.t_1, variables_s, domain)
        erap = ERA5p(onsala.t_0, onsala.t_1, variables_p, domain)

        print (erap.era.time)

        return erap


if __name__ == "__main__":



    variables_p = ['temperature',
                     'u_component_of_wind',
                     'v_component_of_wind',
                     'relative_humidity']

#    variables_s = ["surface_pressure",
                   # "orography",
                   # "skin_temperature",
                   # "2m_temperature",
                   # "10m_u_component_of_wind",
                   # "10m_v_component_of_wind",
                   # "sea_ice_cover",
                   # "land_sea_mask",
                   # "snow_depth"]



    # pressure grid
    p_grid_fine = alt2pres(np.arange(-700, 8000, 125))
    p_grid_coarse = alt2pres(np.arange(8000, 20000, 250))
    p_grid = (np.concatenate([p_grid_fine, p_grid_coarse,
                             np.array([30, 20, 10, 7, 5, 3, 2, 1]) * 100]))

    t0 = datetime(2020, 10, 9, 3, 00, 00 )
    t1 = datetime(2020, 10, 9, 6, 00, 00)

    # location A
    lat = 59.86
    lon = 17.76

    # location B
    #lat = 59.50
    #lon = 18.06

    t2 = t1
    while t0 < t1:

         t0 = t0 + timedelta(minutes = 60)
         loc = locations(lat, lon, t0)
         loc.t_1 = t0 + timedelta(minutes =  60)
         #loc.t_1 = t1
         print (t0.day, t0.hour)

         outpath = os.path.expanduser("~/Documents/arlanda_ERA_data/")

         era_output = run_all_cases(p_grid, loc, loc.t_0, outpath)

         shortnames = era_output.shortname

         grids = dict.fromkeys(shortnames, None)
         for shortname in shortnames:

             grids[shortname] = era_output.interpolate(loc, shortname, method = "nearest")


         with open("arlanda_era5_" + str(t0.day) + str(t0.hour) + "_A.pickle", "wb") as f:
             pickle.dump(grids, f)
             pickle.dump(era_output.era.level.data, f)
         f.close()

         if t0 == t1:
            break
