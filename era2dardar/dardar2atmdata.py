#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 14:18:24 2020

@author: inderpreet
"""

import os
import numpy as np
from era2dardar.DARDAR import DARDARProduct
from era2dardar.atmData import atmdata
from era2dardar.utils.alt2pressure import alt2pres
import typhon.arts.xml as xml
from era2dardar.utils.Z2dbZ import Z2dbZ
from era2dardar.atmData import Z2atmdata



def dardar2atmdata(dardar, cloudsat, erap, eras, p_grid, domain = None):
    """
    This method interpolates different fields to DARDAR grid and saves them
    to be in ARTS xml format.
    
    Only DARDAR filrname is to be provided, 
    the corresponding ERA5 data is selected, downloaded, interpolated

    Parameters
    ----------
    dardarFile : DARDAR class instance

    p_grid     : np.array containing pressure levels for ARTS data [Pa]

    Returns
    -------
    None.

    """

    atm             = atmdata(dardar, erap, eras, p_grid, domain = domain)
    
    time            = atm.t_0.strftime("%Y%m%d%H%M")
    
    vmr_h2o         = atm.vmr_h2o
    vmr_N2          = atm.vmr_N2
    vmr_O2          = atm.vmr_O2
    vmr_O3          = atm.vmr_O3
    lwc             = atm.clwc   
    vmr_field       = np.concatenate([vmr_N2,vmr_O2,
                              vmr_h2o, vmr_O3, lwc], axis = 0)    
    
    
    z_field         = atm.z_field
 
    
    # select reflectivities from cloudsat
    Z2atm           = Z2atmdata(cloudsat, p_grid, z_field)
    
    Z               = Z2atm.Z
    

# atm_fields as a dictionary    
    atm_fields      = {
                         "time"            : time,  
                         "lon_grid"        : atm.lon,
                         "lat_grid"        : atm.lat,
                         "p_grid"          : atm.p_grid,
                         "t_field"         : atm.temperature,
                         "z_field"         : z_field,
                         "lwc"             : atm.clwc,
                         "skt"             : atm.skin_temperature,
                         "t2m"             : atm.t2m,
                         "wind_speed"      : atm.wind_speed,
                         "wind_direction"  : atm.wind_direction,
                         "iwc"             : atm.iwc(),
                         #"z_surface"       : atm.z_surface,
                         "abs_species"     : atm.abs_species,                   
                         "vmr_field"       : vmr_field, 
                         "lsm"             : atm.lsm,
                         "sea_ice_cover"   : atm.sea_ice_cover,
                         "snow_depth"      : atm.snow_depth,
                         "reflectivities"  : Z,
                        # #"N0star"          : atm.N0star,
                         "z_surface"       : atm.z_surface_srtm
                        }

        
    return atm_fields     
        
    
