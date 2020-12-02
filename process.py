#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script matches DARDAR and ERA5 data in time and 
interpolates ERA5 to DARDAR locations

Created on Tue Dec  1 14:43:42 2020

@author: inderpreet
"""
import os
import numpy as np
from datetime import datetime, timedelta
from read_DARDAR_hdf4 import DARDARProduct
from read_ERA import ERA
from dardar2era import dardar2era

import matplotlib.pyplot as plt


filename = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310231507_50677.hdf")

# create DARDAR instance
dardar = DARDARProduct(filename)
  
# time stamp of DARDAR data 
t_0 = dardar.filename_to_date()
t_1 = t_0 + timedelta(minutes = 30)

# create ERA5 product instance based on DARDAR timestamp
ERA = ERA(t_0, t_1)

# regrid ERA5 to DARDAR locations

p_grid = np.arange(1, 1150, 10)
grid_t, grid_z = dardar2era(dardar, ERA, p_grid)




