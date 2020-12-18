#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 12:50:51 2020

@author: inderpreet
"""
import numpy as np

def shift_lon(lon):
    """
    shifts longitudes to the range -[-180, 180] or [0, 360]  to avoid the
    discontinuity arising at 360/0 or 180/-180
    

    Parameters
    ----------
    lon : np.array containing longitude values

    Raises
    ------
    Exception
        if longitude values are smaller than -180 or greater than 360

    Returns
    -------
    np.array containing longitude values in the range avoiding the discontinuity. 

    """
    
    if np.any(lon < -180.0) or np.any(lon > 360.0):
        raise Exception("Longitude values not in valid range")
        
        
    
    lon_360     = lon_180_to_360(lon)
    diff_360    = lon_360.max() - lon_360.min()
    
    lon_180     = lon_360_to_180(lon)
    diff_180    = lon_180.max() - lon_180.min()

    if diff_360 > diff_180:
        
        return lon_180
    else:
        return lon_360
    
    

def lon_180_to_360(lon):
    """
    converts longitudes from range [-180, 180] to [0, 360]

    Parameters
    ----------
    lon : np.array  or scaler
    Returns
    -------
    np.array or scaler in the [0, 360 ] range

    """
    return lon % 360
        

def lon_360_to_180(lon):
    """
    converts longitudes from range [[0, 360] to 

    Parameters
    ----------
    lon : np.array  or scaler
    Returns
    -------
    np.array or scaler in the [-180, 180]  range

    """
    return(lon + 180) % 360 - 180