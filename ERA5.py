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
import numpy as np
from alt2pressure import alt2pres
from interpolator import interpolator


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
             
       
    def expand_lon(self, lon, A):
        """
        Expands longitude of ERA5 by one point to mimic 
        wrapping of data in interpolation
        
        extra longitudnal point 360.0 is added, the corresponding value is copied
        from longitude 0 deg.
        
        Parameters
        ----------
        lon : np.array containing longitude values from 0 - 359.75
        A : np.array containing the values at grid points defined by lon 
        Returns
        -------
        lon : np.array with extended longitudes 
        A : np.array with extra dimensions 
        """
            
        A_repeat = A[:, :, :1]
        A = np.concatenate([A, A_repeat], axis = 2)
    
        lon = np.concatenate((lon, np.arange(0, 0.25, 0.25) + 360.0))
    
        return lon, A     
       
    def interpolate(self, dardar, p_grid = None):
        """
        

        Parameters
        ----------
        dardar : Instance of DARDAR class
        p_grid : if defined, pressure grid for interpolation (hPa) is used
                 otherwise, DARDAR vertical locations are used.
 
        Returns
        -------
        grid_t : np.array of gridded ERA5 data on DARDAR grid

        """
#   get DARDAR locations       
        lon_d       = dardar.get_data('longitude')
        lat_d       = dardar.get_data('latitude')
        height_d    = dardar.get_data('height')
        
        
#   convert longitude from -180 -- 180 to 0--360, if required
        if lon_d.min()  < 0:
            lon_d = lon_d % 360

#   if on pressure levels, read in levels and        
#   add extra pressure level in ERA5 data
        xlevel = 1200 # hPa
        self.add_extra_level(xlevel)
        
        level   = self.era['level'].data 
        level   = np.log(level)  # convert pressure to log            
            
#   get ERA lat/lon/pressure grids and corresponding field    
        lat     = self.era['latitude'].data
        lon     = self.era['longitude'].data        
        field   = self.era[self.shortname].data[0] # 0 for time dimension 
        
#   add one extra dimension to longitude to wrap around during interpolation    
        lon, field = self.expand_lon(lon, field)        
             
#   interpolate ERA5 to DARDAR lat/lon locations    
 
        if p_grid is not None: # if a new p_grid is required  
          
            p = np.log(p_grid) # convert pressure to log range
        else:    # use DARDAR vertical grid

            p_grid = alt2pres(height_d) * 0.01 # convert to hPa
            p = np.log(p_grid)

        pts = []
        for i in range(len(p)):
            point = [[p[i], lat_d[j], lon_d[j]] for j in range(len(lat_d))] 
            pts.append(point)       
                 
        my_interpolating_function = interpolator((level, lat, lon), field)         
        grid_t = my_interpolating_function(pts)
        
        return grid_t
      
        
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
            "ERA5/reanalysis-era5-single-levels/reanalysis-era5-single-levels_"
            + str(year) + str(month) 
            + str(day) + str(hour) + "_"  + parameter + '.nc'])
        
#       load ERA5 data into an xarray           
        self.era = data.open(filename = file[0])        

#       flip latitude to be in ascending order        
        self.era = self.era.sortby('latitude' , ascending = True)     
        
    def expand_lon(self, lon, A):
        """
        Expands longitude of ERA5 by one point to mimic 
        wrapping of data in interpolation
        
        extra longitudnal point 360.0 is added, the corresponding value is copied
        from longitude 0 deg.
        
        Parameters
        ----------
        lon : np.array containing longitude values from 0 - 359.75
        A : np.array containing the values at grid points defined by lon 
        Returns
        -------
        lon : np.array with extended longitudes 
        A : np.array with extra dimensions 
        """
            
        A_repeat = A[:, :1]
        A = np.concatenate([A, A_repeat], axis = 1)
    
        lon = np.concatenate((lon, np.arange(0, 0.25, 0.25) + 360.0))
    
        return lon, A             
         
        
    def interpolate(self, dardar):
        """
        

        Parameters
        ----------
        dardar : Instance of DARDAR class

        Returns
        -------
        grid_t : np.array of gridded surface ERA5 data on DARDAR grid

        """
        lon_d = dardar.get_data('longitude')
        lat_d = dardar.get_data('latitude')
        
        
#   convert longitude from -180 -- 180 to 0--360, if required
        if lon_d.min()  < 0:
            lon_d = lon_d % 360
            
            
#   get ERA lat/lon/pressure grids    
        lat     = self.era['latitude'].data
        lon     = self.era['longitude'].data
        
        field   = self.era[self.shortname].data[0] # 0 for time dimension 
        
#   add one extra dimension to longitude to wrap around during interpolation    
        lon, field = self.expand_lon(lon, field)        

        pts = [[lat_d[j], lon_d[j]] for j in range(len(lat_d))]       
                 
        my_interpolating_function = interpolator((lat, lon), field)         
        grid_t = my_interpolating_function(pts)
        
        return grid_t
        
        
     

