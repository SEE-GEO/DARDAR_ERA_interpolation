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



def dardar2atmdata(dardar, p_grid):
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

    # DARDAR file

#    filename        = dardarfile
    
    # create DARDAR instance
#    dardar          = DARDARProduct(filename, latlims = lat_lims, node = N) 
    
#    p_grid          = alt2pres(np.arange(-100, 25000, 250)) #[Pa] 
    atm             = atmdata(dardar, p_grid)
    vmr_h2o         = atm.vmr_h2o
    vmr_N2          = atm.vmr_N2
    vmr_O2          = atm.vmr_O2
    lwc             = atm.clwc
    vmr_field       = np.concatenate([vmr_N2,vmr_O2,
                              vmr_h2o, lwc], axis = 0)
    time            = atm.t_0.strftime("%Y%m%d%H%M")
    
  
# atm_fields as a dictionary    
    atm_fields      = {
                        "time"            : time,  
                        "lon_grid"        : atm.lon,
                        "lat_grid"        : atm.lat,
                        "p_grid"          : atm.p_grid,
                        "t_field"         : atm.temperature,
                        "z_field"         : atm.z_field,
                        "lwc"             : atm.clwc,
                        "skt"             : atm.skin_temperature,
                        "t2m"             : atm.t2m,
                        "wind_speed"      : atm.wind_speed,
                        "wind_direction"  : atm.wind_direction,
                        "iwc"             : atm.iwc,
                         #Z               : atm.Z,
                        "z_surface"       : atm.z_surface,
                        "abs_species"     : atm.abs_species,                   
                        "vmr_field"       : vmr_field, 
                        "lsm"             : atm.lsm,
                        "sea_ice_cover"   : atm.sea_ice_cover,
                        "snow_depth"      : atm.snow_depth,
                        }
# save  atm_fields to xml files   
    
    # for key in atm_fields.keys():
        
    #     filename = "xmls/" + key + ".xml"
        
    #     xml.save(atm_fields[key], filename)
        
    return atm_fields     
        
    
