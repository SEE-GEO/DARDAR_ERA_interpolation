#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 12:10:13 2020

@author: inderpreet
"""


import os
import numpy as np
from datetime import datetime, timedelta
from era2dardar.DARDAR import DARDARProduct
from era2dardar.ERA5 import ERA5p, ERA5s
from era2dardar.utils.alt2pressure import alt2pres
import typhon.physics.thermodynamics as thermodynamics 
import typhon.physics.atmosphere as atmosphere
import typhon.arts.xml as xml
from typhon.topography import SRTM30
from era2dardar.utils.scale_vmr import scale_vmr
import matplotlib.pyplot as plt

# pressure grid
p_grid = alt2pres(np.arange(-700, 20000, 250))
p_grid = np.concatenate([p_grid, np.array([30, 20, 10, 7, 5, 3, 2, 1]) * 100])
p_grid = p_grid * 0.01



# DARDAR file
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310114257_50670.hdf")
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071190249_47194.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071104824_47189.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/06/07/DARDAR-CLOUD_v2.1.1_2015158155052_48459.hdf")
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2009/06/03/DARDAR-CLOUD_v2.1.1_2009154021055_16481.hdf")
filename = file

# create DARDAR instance
dardar = DARDARProduct(filename, latlims = [-80, 80], node = "A")
  
lon          = dardar.longitude 
lat          = dardar.latitude

lat1         = np.around(lat.min() - 2)
lat2         = np.around(lat.max() + 2)

lon1         = np.around(lon.min() - 2)
lon2         = np.around(lon.max() + 2)

if lon1 < -180:
    lon1 = -180.0
if lon2 > 180:
    lon2 = 180.0    
    
domain  = [lat1, lat2, lon1, lon2]
        
# time stamp of DARDAR data 
t_0 = dardar.t_0
#t_1 = t_0 + timedelta(seconds = float(dardar.time[-1]))
t_1 = t_0 + timedelta(minutes = 30)

# coordinates of DARDAR data
lat_d       = dardar.get_data('latitude')
lon_d       = dardar.get_data('longitude')
lon_d       = lon_d % 360
z_d         = dardar.get_data('height')

# temperature
var         = "specific_humidity"
ERA_t       = ERA5p(t_0, t_1, var, domain = None)
grid_t      = ERA_t.interpolate(dardar, p_grid)


#--------------------------------------------------------------

ERA_tt       = ERA5p(t_0, t_1, var, domain = domain)
t_new       = ERA_tt.interpolate(dardar, p_grid)
 

fig, ax = plt.subplots(1, 1)
ax.plot( lat_d, grid_t[25, :])
ax.plot( lat_d, t_new[25 ,:])

ax.set_yscale('log')
 
#-----------------------------------------------------------
var         = "2m_temperature"
ERA_t       = ERA5s(t_0, t_1, var, domain = None)
grid_t      = ERA_t.interpolate(dardar,)


#--------------------------------------------------------------

ERA_tt       = ERA5s(t_0, t_1, var, domain = domain)
t_new       = ERA_tt.interpolate(dardar)

fig, ax = plt.subplots(1, 1)
ax.plot( grid_t[ :], '.')
ax.plot( t_new[:], '.')



#-------------------------------------------------------------

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
grid_skt      = np.expand_dims(grid_skt, 1)

# 2 m temperature
var           = "2m_temperature"
ERA_t2        = ERA5s(t_0, t_1, var)
grid_t2       = ERA_t2.interpolate(dardar)
grid_t2       = np.expand_dims(grid_t2, 1)


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
wind_speed    = np.expand_dims(wind_speed, 1)
# convert to wind direction
wind_dir      = np.arctan2(grid_u, grid_v) * 180 / np.pi  # degree
wind_dir    = np.expand_dims(wind_dir, 1)


#-------------------------write to xmls---------------------------------
grid_t      = np.expand_dims(grid_t, 2)
grid_z      = np.expand_dims(z_field, 2)
vmr_grid    = np.expand_dims(vmr_grid, 3)

xml.save(p_grid * 100, 'xmls/p_grid.xml')  # Pa
xml.save(lat_d, 'xmls/lat_grid.xml')
xml.save(lon_d, 'xmls/lon_grid.xml')
xml.save(grid_t, 'xmls/t_field.xml')
xml.save(z_field, "xmls/z_field.xml")
xml.save(z_surface, 'xmls/z_surface.xml')
xml.save(grid_skt, 'xmls/skt.xml')
xml.save(wind_speed, 'xmls/wind_speed.xml')
xml.save(wind_dir, 'xmls/wind_direction.xml')
xml.save(abs_species, "xmls/abs_species.xml")
xml.save(grid_t2, 'xmls/t2m.xml')


