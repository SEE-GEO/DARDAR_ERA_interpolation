#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 11:06:11 2020

@author: inderpreet
"""
from scipy.interpolate import RegularGridInterpolator

def interpolate(level, lat, lon, A):
    """
    

    Parameters
    ----------
    level : np.array of dimension N, gives the z coordinate
    lat : np.array of dimension M
    lon : np.array of dimension M
    A : field values to be interpolated from, dimension NxM

    Returns
    -------
    interpolator function

    """
    
    return RegularGridInterpolator((level, lat, lon), A)