#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 12:33:27 2020

@author: inderpreet
"""


import xarray
import os
from datetime import datetime, timedelta

class DARDARProduct():
    
    def __init__(self, filename):
        self.file = filename
        self.dataset = xarray.open_dataset(filename)
        
        
    def get_latitude(self):
        
        return self.dataset["latitude"]
    
    def get_longitude(self):
        
        return self.dataset["longitude"]
    
    def get_height(self):
        
        return self.dataset["height"]
    
    def get_iwc(self):
        
        return self.dataset["iwc"]
    
    def get_temperature(self):
        
        return self.dataset["temperature"]
        
   
    def filename_to_date(self):
        
        filename = os.path.basename(self.file)
        filename = filename.split("_")[2]
        pattern = "%Y%j%H%M%S"
        return datetime.strptime(filename, pattern)
        
    
    def match_era5(self, variable):
        """
        gives the nearest temporal neighbour of DARDAR data

        Parameters
        ----------
        variable : the geophysical parameter to be downloaded
        for example, temperature, 
        Returns
        -------
        era5_filename : string of nearest collocated ERA5 filename
        """

        
        t_dardar = self.filename_to_date()
        
        if t_dardar.minute > 30:
            t_dardar += timedelta(hours = 1)
        
        year  = t_dardar.year
        month = f"{t_dardar.month:02d}"
        day   = f"{t_dardar.day:02d}"
        hour  = f"{t_dardar.hour:02d}"      

            
        era5_filename = (
            "reanalysis-era5-pressure-levels_" + str(year) + str(month) 
            + str(day) + str(hour) + "_"  + variable + '.nc')
            
        return era5_filename    