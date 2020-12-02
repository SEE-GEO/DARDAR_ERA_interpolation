#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 14:44:07 2020

@author: inderpreet
"""
import numpy as np
from expand_lon import expand_lon
from interpolate import interpolate

def dardar2era(dardar, ERA, p_grid):
    """
    interpolates ERA5 data to DARDAR locations
    and the pressure grid is defined in p_grid

    Parameters
    ----------
    dardar :  DARDARProduct instance
    ERA :  ERA5 instance
    p_grid : a pressure grid in hPa, where the values should to interpolated to
    
    Returns
    -------
    grid_t : temperature gridded to DARDAR locations
    grid_z : geopotential gridded to DARDAR locations

    """
    lon_d = dardar.get_data('longitude')
    lat_d = dardar.get_data('latitude')
    height_d = dardar.get_data('height')

    
#   convert longitude from -180-180 to 0-360
    if lon_d.min()  < 0:
        lon_d = lon_d % 360
 
#   add extra pressure level in ERA5 data
    xlevel = 1200
    ERA.add_extra_level('temperature', xlevel)
    ERA.add_extra_level('geopotential', xlevel)
    
#   get ERA lat/lon/pressure grids
    
    lat     = ERA.t.latitude.data
    lon     = ERA.t.longitude.data
    level   = ERA.t.level.data 
    t       = ERA.t.t[0].data
    z       = ERA.z.z[0].data
    
    level   = np.log(level)  # convert pressure to log
    
    
#   add two extra dimension to longitudes to wrap around during interpolation
    
    lon, z = expand_lon(ERA.z.longitude.data, z )
    lon, t = expand_lon(ERA.t.longitude.data, t )
            
    #my_interpolating_function = RegularGridInterpolator((level, lat, lon), A)
    
    p_grid = np.arange(1, 1150, 10)
    points = []
    
#   interpolate ERA5 to DARDAR lat/lon locations
    
    for i in range(len(p_grid)):
        p = np.log(p_grid[i]) # convert pressure to log range
        pts = [[p, lat_d[j], lon_d[j]] for j in range(len(lat_d))]       
        points.append(pts)
             
    my_interpolating_function = interpolate(level, lat, lon, t)         
    grid_t = my_interpolating_function(points)
    
    my_interpolating_function = interpolate(level, lat, lon, z)         
    grid_z = my_interpolating_function(points)
    
    return grid_t, grid_z