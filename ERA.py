#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 08:58:25 2020

Interface to download and load ERA5 data nearest to DARDAR time
Allows to add extra pressure layer below 1000 hPa

The download and loading of ERA5 data uses pansat Product class "ERA5Product"

@author: inderpreet
"""

from pansat.products.reanalysis.era5 import ERA5Product
import xarray
from alt2pressure import pres2alt
from scipy.constants import g

class ERA():
    """
    Downloads the required ERA5 data as per DARDAR timestamp and load the 
    required fields
    """
    
    
    def __init__(self, t_0, t_1):
        """
        download and read ERA5 data between times t_0 and t_1

        Parameters
        ----------
        t_0 : datetime.datetime object, start time
        t_1 : datetime.datetime object, end time

        Returns
        -------
        None.

        """

        temp = ERA5Product('hourly','pressure', ['temperature'])
        z    = ERA5Product('hourly','pressure', ['geopotential']) 
  
        file_temp = temp.download(t_0, t_1)
        file_z    = z.download(t_0, t_1)
        
        year  = t_1.year
        month = f"{t_1.month:02d}"
        day   = f"{t_1.day:02d}"
        hour  = f"{t_1.hour:02d}"      

#        file_temp = "ERA5/reanalysis-era5-pressure-levels/reanalysis-era5-pressure-levels_2020030112_temperature.nc"
#        file_z    = "ERA5/reanalysis-era5-pressure-levels/reanalysis-era5-pressure-levels_2020030112_geopotential.nc"
        file_temp = ([
            "ERA5/reanalysis-era5-pressure-levels/reanalysis-era5-pressure-levels_"
            + str(year) + str(month) 
            + str(day) + str(hour) + "_"  + "temperature" + '.nc'])
        
        file_z = ([
            "ERA5/reanalysis-era5-pressure-levels/reanalysis-era5-pressure-levels_"
            + str(year) + str(month) 
            + str(day) + str(hour) + "_"  + "geopotential" + '.nc'])
            
        print (file_temp, file_z)
        self.t = temp.open(filename = file_temp[0])        
        self.z    = z.open(filename = file_z[0])

# flip latitude to be in ascending order        
        self.t = self.t.sortby('latitude' , ascending = True) 
        self.z = self.z.sortby('latitude' , ascending = True)
        
    
    def add_extra_level(self, variable, xlevel):
         """
         adds extra pressure level at xlevel hPa to allow interpolation to 
         levels lower than 1000 hPa.
         Temperature at xlevel hPa is equal to 1000 hPa
         geopotential is converted using hydrostatic equation
         
         Parameters
        ----------
        variable = geophysical parameter for which extra pressure layer 
        is added 
        xlevel = extra pressure level in hPa
         Returns
         -------
         None.
         
         """
         
         if variable not in ["geopotential", "temperature"]:
             raise Exception("variable should be one of [geopotential,temperature]")
             
         if variable == "geopotential":
          # geopotential                    
             A          = self.z.z[:, -1, :, :].to_dataset() # copy lowest pressure level
             A["level"] = xlevel
             self.z     = (xarray.concat([self.z, A], dim="level"))
             
             # convert pressure to geopotential
             self.z.z[0, -1, :, :] = pres2alt(xlevel * 100) * g
         
         else:   
             # temperature
             A          = self.t.t[:, -1, :, :].to_dataset() # copy lowest pressure level
             A["level"] = xlevel
             self.t     = (xarray.concat([self.t, A], dim="level"))
             
         
         
     

