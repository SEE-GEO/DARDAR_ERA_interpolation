#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 11:27:59 2020

@author: inderpreet
"""
from typhon import constants

def ellipsoidmodels(model):
    """
    ccdcs

    Parameters
    ----------
    model : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    ellipsoid : TYPE
        DESCRIPTION.

    """
    
    
    if model == 'sphericalearth':
        ellipsoid = [ constants.earth_radius, 0 ]

    elif model == 'wgs84':
         ellipsoid = [ 6378137, 0.081819190842621 ]
    
    else:
        raise Exception("Unknown model, only spherical earth and wgs84 are supported")
        
    return ellipsoid    
     

    