#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 2 08:58:25 2020

Interface to download and load ERA5 data nearest to DARDAR time
Adjusts ERA5 data to increasing latitudes
and also
Allows to add extra pressure layer below 1000 hPa with 
function "add_extra_level" 

The download and loading of ERA5 data uses pansat Product class "ERA5Product"
See README to install pansat.

@author: inderpreet
"""
from era2dardar.ERA5_parameters import parameters
from pansat.products.reanalysis.era5 import ERA5Hourly
import xarray
from era2dardar.utils.alt2pressure import pres2alt, alt2pres
from scipy.constants import g
import numpy as np
from era2dardar.utils.interpolator import interpolator



class ERA5():
    """
    Downloads the required ERA5 data as per DARDAR timestamp and load the 
    required fields
    """
    
    
    def __init__(self, t_0, t_1,  variables, domain = None):
        """
        download and read ERA5 data between times t_0 and t_1

        Parameters
        ----------
        t_0 : datetime.datetime object, start time
        t_1 : datetime.datetime object, end time
        variables : list containing the longnames of ERA5 fields 
        to be downloaded, see  ERA5_parameters.parameters for listed parameters
        """
        
        # if parameter not in parameters:
        #     raise Exception(
        #         "parameter currently not listed check spelling or add to parameters dictionary", 
        #                     parameters)
        self.longname = []
        self.shortname = []
        for parameter in variables: 
            self.longname.append(parameter)
            self.shortname.append(parameters[parameter])  
            self.t_0 = t_0
            self.t_1 = t_1


        
class ERA5p(ERA5):
    """
    class to download and load ERA5 surface variables
    inherits from class ERA5
    """
    
    def __init__(self, t_0, t_1, variables, domain = None):
        
        super().__init__(t_0, t_1, variables, domain = None)

        self.domain = domain

        # create ERA5product instance
        data = ERA5Hourly('pressure', self.longname, domain = self.domain)
  
        # download data with matching time stamp
        file = data.download(self.t_0, self.t_1)
        
        # load ERA5 data into an xarray           
        self.era = data.open(filename = file[0])   

        # flip latitude to be in ascending order        
        self.era = self.era.sortby('latitude' , ascending = True) 
             
        # add extra pressure level in ERA5 data
        xlevel = 1150 # hPa
        self.add_extra_level(xlevel)
        

        
        
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
         # extract dataset at last pressure level
         Ex             =  self.era[dict(level = [-1])]
         
         # geopotential at 1150 hPa is using simple conversion
         if "geopotential" in self.longname:
             Ex["z"]        = pres2alt(xlevel * 100) * g
             
         Ex["level"]    = np.array([xlevel])
         self.era   = (xarray.concat([self.era, Ex], dim="level"))        
         
       
    def expand_lon(self, lon, A):
        """
        Expands longitude of ERA5 by one point to mimic 
        wrapping of data in interpolation
        
        extra longitudnal point 180.25 and -180.25 are added, the corresponding value is copied
        from longitude -180.0 and 180.0 deg respectively
        
        Parameters
        ----------
        lon : np.array containing longitude values from -180 to 180
        A : np.array containing the values at grid points defined by lon 
        Returns
        -------
        lon : np.array with extended longitudes 
        A : np.array with extra dimensions 
        """

        A_start = A[:, :, -1:]
        A_end   = A[:, :, :1]
        A = np.concatenate([A_start, A, A_end], axis = 2)
 
        lon = np.concatenate(([lon.min() - 0.25], lon))    
        lon = np.concatenate((lon, [0.25 + lon.max()]))

        return lon, A     
       
    def interpolate(self, other, shortname, p_grid = None, method = "linear"):
        """
        

        Parameters
        ----------
        other : Instance of DARDAR/locations class
        p_grid : if defined, pressure grid for interpolation (hPa) is used
                 otherwise, DARDAR vertical locations are used.
 
        method : "nearest", "linear"; default is "linear"
        Returns
        -------
        grid_t : np.array of gridded ERA5 data on DARDAR grid

        """
        #   get DARDAR locations       
        lon_d       = other.longitude
        lat_d       = other.latitude
        


        
        level   = self.era['level'].data 
        level   = np.log(level)  # convert pressure to log            
            
        #   get ERA lat/lon/pressure grids and corresponding field    
        lat     = self.era['latitude'].data
        lon     = self.era['longitude'].data     
        field   = self.era[shortname].data[0] # 0 for time dimension 

        #   add one extra dimension to longitude to wrap around during interpolation
        if lon.min() == -180.0:
            lon, field = self.expand_lon(lon, field)      
        
        #convert to 0 to 360, this is needed when global data is downloaded, domain = None
        if lon.max() > 180.5:
            lon_d  = lon_d % 360.0
            
        #interpolate ERA5 to DARDAR lat/lon locations 
 
        if p_grid is not None: # if a new p_grid is required           
            p = np.log(p_grid) # convert pressure to log range

        else:    # use ERA5 vertical grid
            p = level

        pts = []
        for i in range(len(p)):
            point = [[p[i], lat_d[j], lon_d[j]] for j in range(len(lat_d))] 
            pts.append(point)       
                 
        my_interpolating_function = (interpolator((level, lat, lon),
                                                  field, method = method))        
        grid_t = my_interpolating_function(pts)
        
        return grid_t
      
        
class ERA5s(ERA5):
    """
    class to download and load ERA5 surface variables
    inherits from class ERA5
    """
    
    def __init__(self, t_0, t_1, parameter, domain= None):
        
        super().__init__(t_0, t_1, parameter, domain = None)
        
        self.domain = domain

        #  create ERA5product instance
        data = ERA5Hourly('surface', self.longname, domain = self.domain)
  
        # download data with matching time stamp
        file = data.download(self.t_0, self.t_1)
        
         
        self.era = data.open(filename = file[0])        

        #  flip latitude to be in ascending order        
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
        
        A_start = A[ :, -1:]
        A_end   = A[ :,  :1]
        A = np.concatenate([A_start, A, A_end], axis = 1)
 
        lon = np.concatenate(([lon.min() - 0.25], lon))    
        lon = np.concatenate((lon, [0.25 + lon.max()]))
    
        return lon, A             
         
        
    def interpolate(self,  other, shortname, method = "linear"):
        """
        

        Parameters
        ----------
        other : Instance of DARDAR/locations class

        Returns
        -------
        grid_t : np.array of gridded surface ERA5 data on DARDAR grid
        method : "linear", "nearest", default is "linear"
        
        """
        lon_d = other.longitude
        lat_d = other.latitude
        
        

            
            
        #   get ERA lat/lon/pressure grids    
        lat     = self.era['latitude'].data
        lon     = self.era['longitude'].data
        
        field   = self.era[shortname].data[0] # 0 for time dimension 
        

        #   add one extra dimension to longitude to wrap around during interpolation  
        if lon.min() == -180.0:
            lon, field = self.expand_lon(lon, field)        

        #   convert longitude from -180 -- 180 to 0--360, if required
        #   in some cases ERA5 longitudes could be 0 to 360 format            
        if lon.max() > 180.5:
            lon_d  = lon_d % 360.0
            
        pts = [[lat_d[j], lon_d[j]] for j in range(len(lat_d))]       
        
        my_interpolating_function = (interpolator((lat, lon), field, 
                                                  method= method))        
        grid_t = my_interpolating_function(pts)
        
        return grid_t
        
        
     

