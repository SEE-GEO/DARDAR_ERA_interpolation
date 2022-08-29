#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 11:04:58 2021

This script is used to generate the csv files containing refractive index data
provided to Ericsson. This script does not download the ERA5 data, but uses
instead the bulk files (containing data for one month) I downloaded seperately.

@author: inderpreet
"""
import os
import xarray as xr
import numpy as np
import pandas as pd
import glob
from typhon.physics.thermodynamics import e_eq_water_mk
from era2dardar.utils.refractive_index import refractive_index
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from era2dardar.utils.alt2pressure import pres2alt
plt.rcParams.update({'font.size': 16})

# %%
# location A
key = "B"
if key == "A":
    lat = 59.850
    lon = 17.767
if key == "B":
    lat = 59.511
    lon = 18.067

inpath = "/home/inderpreet/data/Ericsson_data/"

files = glob.glob(os.path.join(inpath, "*.nc"))

for file in files:

    filename = os.path.join(inpath, file)

    era5 = xr.open_dataset(filename)

    era5 = era5.interp(latitude=lat, longitude=lon, method="nearest")

    T_era = era5["t"].values
    r_era = era5["r"].values
    p_era = era5["level"].values * 100
    p = np.repeat(p_era.reshape(-1, 1), era5.time.size, axis=1).T
    e_w_era = r_era * e_eq_water_mk(T_era)/100

    n_era = refractive_index(p, e_w_era, T_era)
    times = era5.time.values
    n_era_ppm = (n_era - 1) * 1000000

    era_h = pres2alt(p_era)

    for i in range(times.size):
        print(times[i])
        loc_data = {}

        loc_data["height"] = era_h[:-1]
        loc_data["pressure"] = p_era[:-1]
        loc_data["relative_humidity"] = r_era[i, :-1]
        loc_data["temperature"] = T_era[i, :-1]

        loc_data["refractive_index"] = n_era_ppm[i, :-1]

        df = pd.DataFrame.from_dict(loc_data)

        time = datetime.utcfromtimestamp(times[i].tolist()/1e9)
        local_time = time - timedelta(minutes=60)
        outfile = "loc_" + key + "_" + local_time.strftime("%Y%m%d:%H")
        print(outfile)
        df.to_csv("refractive_index_data_ERA5_2020/" + outfile + ".csv")
