#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 14:03:56 2020

@author: inderpreet
"""
import numpy as np

def expand_lon(lon, A):
    """
    Expands longitude of ERA5 by one point to mimic 
    wrapping of data in interpolation
    
    extra longitudnal point 360.0 is added, the corresponding value is copied
    from longide 0 deg.
    
    Parameters
    ----------
    lon : np.array containing longitude values from 0 - 359.75
    A : np.array containing the values at grid points defined by lon 
    Returns
    -------
    lon : np.array with extended longitudes 
    A : np.array with extra dimensions repeating data 
    """
    
    A_repeat = A[:, :, :1]
    A = np.concatenate([A, A_repeat], axis = 2)

    lon = np.concatenate((lon, np.arange(0, 0.25, 0.25) + 360.0))

    return lon, A