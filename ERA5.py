#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 08:58:25 2020

Interface to download and load ERA5 data nearest to DARDAR time
Adjusts ERA5 data to increasing latitudes
and also
Allows to add extra pressure layer below 1000 hPa with 
function "add_extra_level" 

The download and loading of ERA5 data uses pansat Product class "ERA5Product"
See README to install pansat.

@author: inderpreet
"""
from ERA5_parameters import parameters
from pansat.products.reanalysis.era5 import ERA5Product
import xarray
from alt2pressure import pres2alt
from scipy.constants import g

class ERA5():
    """
    Downloads the required ERA5 data as per DARDAR timestamp and load the 
    required fields
    """
    
    
    def __init__(self, t_0, t_1, parameter):
        """
        download and read ERA5 data between times t_0 and t_1

        Parameters
        ----------
        t_0 : datetime.datetime object, start time
        t_1 : datetime.datetime object, end time
        parameter : string containing the longname of ERA5 field 
        to be downloaded, see  ERA5_parameters.parameters for listed parameters
        """
        
        if parameter not in parameters:
            raise Exception(
                "parameter currently not listed check spelling or add to parameters dictionary", 
                            parameters)
        
        
        self.longname  = parameter
        self.shortname = parameters[parameter]       


class ERA5p(ERA5):
    """
    class to download and load ERA5 surface variables
    inherits from class ERA5
    """
    
    def __init__(self, t_0, t_1, parameter):
        super().__init__(t_0, t_1, parameter)
        
#       create ERA5product instance
        data = ERA5Product('hourly', 'pressure', [parameter])
  
#       download data with matching time stamp
        file = data.download(t_0, t_1)
        
        year  = t_1.year
        month = f"{t_1.month:02d}"
        day   = f"{t_1.day:02d}"
        hour  = f"{t_1.hour:02d}"      

        file = ([
            "ERA5/reanalysis-era5-pressure-levels/reanalysis-era5-pressure-levels_"
            + str(year) + str(month) 
            + str(day) + str(hour) + "_"  + parameter + '.nc'])
        
#       load ERA5 data into an xarray           
        self.era = data.open(filename = file[0])        

#       flip latitude to be in ascending order        
        self.era = self.era.sortby('latitude' , ascending = True) 
    
    def add_extra_level(self,  xlevel):
         """
         adds extra pressure level at xlevel hPa to allow interpolation to 
         levels lower than 1000 hPa for "pressure" ERA5 variables
         Field values at xlevel hPa are equal to 1000 hPa
         however,
         Geopotential is converted using hydrostatic equation
         
         Parameters
        ----------
        variable = geophysical parameter for which extra pressure layer 
        is added 
        xlevel = extra pressure level in hPa
         Returns
         -------
         None.
         
         """         
             
         if self.longname == "geopotential":
          # geopotential                    
             A          = self.era[self.shortname][:, -1, :, :].to_dataset() # copy lowest pressure level
             A["level"] = xlevel
             self.era   = (xarray.concat([self.era, A], dim="level"))
             
             # convert pressure to geopotential
             self.era[self.shortname][0, -1, :, :] = pres2alt(xlevel * 100) * g
         
         else:   
             A          = self.era[self.shortname][:, -1, :, :].to_dataset() # copy lowest pressure level
             A["level"] = xlevel
             self.era   = (xarray.concat([self.era, A], dim="level"))
     
        
class ERA5s(ERA5):
    """
    class to download and load ERA5 surface variables
    inherits from class ERA5
    """
    
    def __init__(self, t_0, t_1, parameter):
        
        super().__init__(t_0, t_1, parameter)
        
#       create ERA5product instance
        data = ERA5Product('hourly', 'surface', [parameter])
  
#       download data with matching time stamp
        file = data.download(t_0, t_1)
        
        year  = t_1.year
        month = f"{t_1.month:02d}"
        day   = f"{t_1.day:02d}"
        hour  = f"{t_1.hour:02d}"      

        file = ([
            "ERA5/reanalysis-era5-pressure-levels/reanalysis-era5-pressure-levels_"
            + str(year) + str(month) 
            + str(day) + str(hour) + "_"  + parameter + '.nc'])
        
#       load ERA5 data into an xarray           
        self.era = data.open(filename = file[0])        

#       flip latitude to be in ascending order        
        self.era = self.era.sortby('latitude' , ascending = True)          
         
     

