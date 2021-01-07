#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 11:58:52 2020
  Earth's gravitational constant
  Returns g as a function of latitude and altitude
  
  follows pos2g from atmlab
@author: inderpreet
"""


import numpy as np
# %
# % FORMAT   g = pos2g( lat [, z] )
# %        
# % OUT   g     Gravitational constant
# % IN    lat   Latitude(s).
# % OPT   z     Altitude(s). Default 0 m.

# % 2006-12-06   Created by Patrick Eriksson.


def pos2g( lat, z = 0):
    """
    

    Parameters
    ----------
    lat : latitudes
    z : altitude, [m] The default is 0.

    Raises
    ------
    Exception if latitude values are outside [-90. 90]
    or altitudes are outside [-1, 1000]km

    Returns
    -------
    None.

    """
   
  
    if np.any( lat<-90 )  or  np.any( lat > 90 ):
        raise Exception( 'Only latitudes inside [-90,90] are allowed.')

    if np.any( z<-1e3 )  or  np.any( z>1000e3 ):
      raise Exception( 'Only altitudes inside [-1,1000] km are allowed.' )
      
    lat = lat * np.pi/180.0
    
     
 
# Expression found on web page of UK's National Physical Laboratory
    
    g = (9.780327 * ( 1 + 0.0053024 * np.sin(lat)**2 
                     - 0.0000058*np.sin(2*lat)**2) - 3.086e-6*z)
#Move to apparent gravity, i.e. include effect of the centrifugal force. See:
# A first course in Atmospheric Thermodynamics by G. Petty (page 89)
# As well as https://glossary.ametsoc.org/wiki/Apparent_gravity
#
    g = g - ((7.29e-5)**2 * 6378e3)*np.cos(lat)**2
                                       
                                                                       
    return g                                                                           