#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 11:34:21 2020

Geocentric radius of a reference ellipsoid

follows the ellipsoidradii function from atmlab

@author: inderpreet
"""
import numpy as np

def ellipsoidradii(ellipsoid, lat):
    """
    

    Parameters
    ----------
    ellipsoid : TYPE
        DESCRIPTION.
    lat : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """

    if len(ellipsoid) < 2:
        raise Exception( "The arg *ellipsoid* must be a vector of length 2" )
    
# Spherical case (all radii the same)
    if ellipsoid[1] == 0:
      r = np.tile( ellipsoid[0], lat.size);
     
    else:
# A re-arranged form of the follwing expression is used:
#       r = ab/sqrt(b^2*cos(lat)^2+a^2*sin(lat)^2)
      
      c = 1 - ellipsoid[1] **2;
      b = ellipsoid[0] * np.sqrt( c );

      lat = lat * np.pi/180.0 # convert to radians
      r = b / np.sqrt( c * np.cos(lat)**2 + np.sin(lat)**2 )
      
      return r
    
