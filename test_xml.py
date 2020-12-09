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
import typhon.physics.thermodynamics as thermodynamics 
import typhon.physics.atmosphere as atmosphere
import typhon.arts.xml as xml
from typhon.topography import SRTM30
from scale_vmr import scale_vmr

p_grid = alt2pres(np.arange(-100, 25000, 250)) * 0.01



# DARDAR file
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310114257_50670.hdf")
filename = file

# create DARDAR instance
dardar = DARDARProduct(filename, latlims = [20, 25])
  
# time stamp of DARDAR data 
t_0 = dardar.filename2date()
t_1 = t_0 + timedelta(minutes = 30)

# coordinates of DARDAR data
lat_d       = dardar.get_data('latitude')
lon_d       = dardar.get_data('longitude')
z_d         = dardar.get_data('height')

# temperature
var         = "temperature"
ERA_t       = ERA5p(t_0, t_1, var)
grid_t      = ERA_t.interpolate(dardar, p_grid)


grid_p      = np.tile(p_grid, (grid_t.shape[1], 1))
grid_p      = grid_p.T 

# z_field (geometrical altitudes fulfilling hydrostatic equation)

z_field     = [[atmosphere.pressure2height(p_grid * 100, grid_t[:, i])] for i in range(grid_t.shape[1])]
z_field     = np.concatenate(z_field).T

# volume mixing ratio
var         = "specific_cloud_liquid_water_content"
ERA_lwc     = ERA5p(t_0, t_1, var)
grid_lwc    = ERA_lwc.interpolate(dardar, p_grid)
rho         = thermodynamics.density(grid_p * 100, grid_t) # density of air
grid_lwc    = grid_lwc * rho  # convert lwc units to mass concentration


var         = "specific_humidity"
ERA_q       = ERA5p(t_0, t_1, var)
grid_q      = ERA_q.interpolate(dardar, p_grid)
grid_q2vmr  = thermodynamics.specific_humidity2vmr(grid_q)

grid_N2     = np.ones(grid_lwc.shape) * 0.781 # vmr N2
grid_O2     = np.ones(grid_lwc.shape) * 0.209 # vmr O2

grid_N2     = scale_vmr(grid_N2, grid_q2vmr)
grid_O2     = scale_vmr(grid_O2, grid_q2vmr)

vmr_grid    = np.stack([grid_N2,
                      grid_O2,
                      grid_q2vmr,
                      grid_lwc])

abs_species   = ["N2","O2","H2O","LWC"]

#  z-surface, from typhon.topography
z_surface     = SRTM30.interpolate(lat_d, lon_d)


# skin temperature
var           = "skin_temperature"
ERA_skt       = ERA5s(t_0, t_1, var)
grid_skt      = ERA_skt.interpolate(dardar)

# 2 m temperature
var           = "2m_temperature"
ERA_t2        = ERA5s(t_0, t_1, var)
grid_t2       = ERA_t2.interpolate(dardar)


# u 10m 
var           = "10m_u_component_of_wind"
ERA_u         = ERA5s(t_0, t_1, var)
grid_u        = ERA_u.interpolate(dardar)

# v 10 m 
var           = "10m_v_component_of_wind"
ERA_v         = ERA5s(t_0, t_1, var)
grid_v        = ERA_v.interpolate(dardar)


# convert to wind speed
wind_speed    = np.sqrt(grid_u**2 + grid_v**2)

# convert to wind direction
wind_dir      = np.arctan2(grid_u, grid_v) * 180 / np.pi  # degree



#-------------------------write to xmls---------------------------------

xml.save(p_grid * 100, 'p_grid.xml')  # Pa
xml.save(lat_d, 'lat_grid.xml')
xml.save(lon_d, 'lon_grid.xml')
xml.save(grid_t, 't_field.xml')
xml.save(z_field, "z_field.xml")
xml.save(z_surface, 'z_surface.xml')
xml.save(grid_skt, 'surface_skin_temperature.xml')
xml.save(wind_speed, 'wind_speed.xml')
xml.save(wind_dir, 'wind_direction.xml')
xml.save(vmr_grid, 'vmr_field.xml')
xml.save(abs_species, "abs_species.xml")


