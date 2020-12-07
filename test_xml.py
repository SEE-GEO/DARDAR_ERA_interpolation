#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 12:10:13 2020

@author: inderpreet
"""


import os
import numpy as np
from datetime import datetime, timedelta
from DARDAR import DARDARProduct
from ERA5 import ERA5p, ERA5s
from alt2pressure import alt2pres
from typhon.physics.thermodynamics import specific_humidity2vmr, density
from typhon.physics.atmosphere import relative_humidity2vmr
import typhon.arts.xml as xml


p_grid = alt2pres(np.arange(-100, 25000, 250)) * 0.01
print(p_grid)


# DARDAR file
file = "~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310231507_50677.hdf"
filename = os.path.expanduser(file)

# create DARDAR instance
dardar = DARDARProduct(filename, latlims = [20, 25])
  
# time stamp of DARDAR data 
t_0 = dardar.filename2date()
t_1 = t_0 + timedelta(minutes = 30)

# temperature
var         = "temperature"
ERA_t       = ERA5p(t_0, t_1, var)
grid_t      = ERA_t.interpolate(dardar, p_grid)


grid_p      = np.tile(p_grid, (grid_t.shape[1], 1))


var         = "specific_cloud_liquid_water_content"
ERA_lwc     = ERA5p(t_0, t_1, var)
grid_lwc    = ERA_lwc.interpolate(dardar, p_grid)

rho         = density(grid_p.T, grid_t) # density of air

grid_lwc    = grid_lwc * rho  # convert lwc units to mass concentration



var         = "specific_humidity"
ERA_q       = ERA5p(t_0, t_1, var)
grid_q      = ERA_q.interpolate(dardar, p_grid)
grid_q2vmr  = specific_humidity2vmr(grid_q)


grid_N2     = np.ones(grid_lwc.shape) * 0.781 # vmr N2
grid_O2     = np.ones(grid_lwc.shape) * 0.209 # vmr O2

grid_data = np.stack([grid_N2,
                      grid_O2,
                      grid_q2vmr,
                      grid_lwc])
xml.save(grid_data, 'vmr_field.xml')


