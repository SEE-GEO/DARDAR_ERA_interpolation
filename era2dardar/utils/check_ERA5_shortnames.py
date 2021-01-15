#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 14:17:42 2021

@author: inderpreet
"""
from pansat.products.reanalysis.era5 import ERA5Product
from datetime import datetime, timedelta

def  check_ERA5p_shortnames(parameter):  
    """
    simple function to check the shortnames of ERA5 pressure variables

    Parameters
    ----------
    parameter : string containing longname of pressure level variable

    Returns
    -------
    xarray describing contents of the ERA5 field

    """    
    
    try:
        data = ERA5Product('hourly', 'pressure', [parameter])        
    except:
        print ("maybe single level")  
    try:
        data = ERA5Product('hourly', 'surface', [parameter])  
    except:    
        raise ValueError("check name of the parameter")
        
  
#       download data with matching time stamp
    t_0 = datetime(2019, 1, 1, 1)
    t_1 = t_0 + timedelta(minutes = 30)
    file = data.download(t_0, t_1)
    era = data.open(filename = file[0])   
    
    return era

