#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 14:44:07 2020

Regrids ERA5 data to DARDAR grid

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
    p_grid : np.array, pressure grid in hPa, where the values should to interpolated to
    parameter : string containing the name of ERA5 field to be interpolated
    Returns
    -------
    grid_t : temperature gridded to DARDAR locations
    grid_z : geopotential gridded to DARDAR locations

    """  
#   read in DARDAR locations        
    lon_d = dardar.get_data('longitude')
    lat_d = dardar.get_data('latitude')

    
#   convert longitude from -180 -- 180 to 0--360
    if lon_d.min()  < 0:
        lon_d = lon_d % 360
 
#   add extra pressure level in ERA5 data
    xlevel = 1200 # hPa
    ERA.add_extra_level(xlevel)
    
#   get ERA lat/lon/pressure grids    
    lat     = ERA.era['latitude'].data
    lon     = ERA.era['longitude'].data
    level   = ERA.era['level'].data 
    field   = ERA.era[ERA.shortname].data[0]

    
    level   = np.log(level)  # convert pressure to log
    
    
#   add two extra dimension to longitudes to wrap around during interpolation    
    lon, field = expand_lon(lon, field)
        
#   interpolate ERA5 to DARDAR lat/lon locations    
    points = []
    
    for i in range(len(p_grid)):
        p = np.log(p_grid[i]) # convert pressure to log range
        pts = [[p, lat_d[j], lon_d[j]] for j in range(len(lat_d))]       
        points.append(pts)
             
    my_interpolating_function = interpolate(level, lat, lon, field)         
    grid_t = my_interpolating_function(points)
        
    return grid_t