#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:45:57 2021

a small class to define lat/lon locations for one or more locations e.g. Onsala

@author: inderpreet
"""

from datetime import datetime, timedelta
import numpy as np

class locations():
    
    def __init__(self, lat, lon, time):
        """
        

        Parameters
        ----------
        lat : list or scalar, the latitude for all locations 
        lon : list or scalar, the longitude for all locations
        time : datetime object

        Raises
        ------
        Exception
            DESCRIPTION.is size of lat and lon are not equal
        Returns
        -------
        None.

        """
        
      
        if np.isscalar(lat):
            lat = [lat]
            lon = [lon]
        
        if len(lat) != len(lon):
            raise Exception("size of lat and lon should be equal")
        
        
        self.latitude  = np.array(lat)
        self.longitude = np.array(lon)
            
            
            
        self.t_0 = time
        self.t_1 = time + timedelta(minutes = 30)