#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 09:18:07 2020

class to convert ERA5 data to DARDAR grid

the scripts downloads the matching ERA5 data in time 
and 
functions for different variables are provided.
The returned outputs for each are in ARTS grid format and can be saved to xml
using typhon xml functions

@author: inderpreet

"""
import os
import numpy as np
from datetime import datetime, timedelta
from era2dardar.ERA5 import ERA5p, ERA5s
from era2dardar.utils.alt2pressure import alt2pres
import typhon.physics.thermodynamics as thermodynamics 
import typhon.physics.atmosphere as atmosphere
import typhon.arts.xml as xml
from typhon.topography import SRTM30
from era2dardar.utils.scale_vmr import scale_vmr
from scipy import interpolate

class atmdata():
    """
    atmdata Class which interpolates ERA5 data to DARDAR/Cloudsat grids
    inputs are dardar class object and the pressure grid over which values are
    to be interpolated
    """
    

    def __init__(self, dardar, p_grid = None):
        """
        
        Parameters
        ----------
        dardar : a DARDAR class instance
        p_grid : np.array, the pressure grid over which ERA5 
        is to be interpolated. Units are in [Pa]
        

        Returns
        -------
        None.

        """        
        
        # time stamp of DARDAR data
        self.dardar  = dardar
        self.time    = dardar.filename2date()
        self.p_grid  = p_grid
        
        self.lat     = dardar.get_data("latitude")

    @property    
    def t_0(self):
        """
        the start time of DARDAR pass

        Returns
        -------
        a datetime object containing the time of DARDAR pass
        """
        return self.time
    
    @property
    def t_1(self):
        """
        the end time of DARDAR pass

        Returns
        -------
        a datetime object containing the end time of DARDAR pass
        """
        return self.time + timedelta(minutes = 30)
    
        
    @property    
    def temperature(self):
        """
        interpolated ERA5 temperature fields to DARDAR grid and pressure grid
        defined in self.p_grid

        Returns
        -------
        grid_t : np.array containing the interpolated values
        dimensions [p, lat, lon]

        """

        p           = self.p_grid        
        var         = "temperature"
        ERA_t       = ERA5p(self.t_0, self.t_1, var)
        grid_t      = ERA_t.interpolate(self.dardar, p)
        grid_t      = np.expand_dims(grid_t, 2)
        
        return grid_t
    
    @property    
    def clwc(self):
        """
        interpolated ERA5 CLWC fields to DARDAR grid and pressure grid
        defined in self.p_grid

        Returns
        -------
        grid_lwc : np.array containing the interpolated values
        dimensions [1, p, lat, lon]

        """
        
        p           = self.p_grid   
            
        var         = "specific_cloud_liquid_water_content"
        ERA_lwc     = ERA5p(self.t_0, self.t_1, var)
        grid_lwc    = ERA_lwc.interpolate(self.dardar, p)
        
        grid_p      = np.tile(p, (grid_lwc.shape[1], 1))
        grid_p      = grid_p.T 
        A           = self.temperature
        rho         = thermodynamics.density(grid_p * 100, A) #air density
        grid_lwc    = grid_lwc * rho  #convert lwc units to mass concentration
        
        grid_lwc    = np.expand_dims(grid_lwc, axis = (0, 3))
        
        return grid_lwc
    
    @property
    def abs_species(self):  
        """
        The absorption species are simply stored as a list of strings

        Returns
        -------
        String containing the names of absorption species

        """        

        abs_species   = ["N2","O2","H2O","LWC"]
        
        return abs_species

    @property
    def z_surface(self):
        """
        z_surface fields interpolated to DARDAR grid.
        The values are interpolated using topography module in typhon
        uses SRTM30 near-global digital elevation model (DEM) 

        Returns
        -------
        z_surface : np.array containing the interpolated values
        dimensions [lat, lon]

        """
        
        #  z-surface, from typhon.topography
        lat_d         = self.dardar.get_data("latitude")
        lon_d         = self.dardar.get_data("longitude")
        z_surface     = SRTM30.interpolate(lat_d, lon_d)
        
        z_surface     = np.expand_dims(z_surface, axis = 1)
        
        return z_surface

    @property
    def skin_temperature(self):
        """
        ERA5 skin temperature fields interpolated to DARDAR grid.
       
        Returns
        -------
        grid_skt : np.array containing the interpolated values
        dimensions [lat, lon]

        """
        var           = "skin_temperature"
        ERA_skt       = ERA5s(self.t_0, self.t_1, var)
        grid_skt      = ERA_skt.interpolate(self.dardar)
        grid_skt      = np.expand_dims(grid_skt, axis = 1)
        
        return grid_skt

    @property
    def t2m(self):
        """
        ERA5 2m temperature fields interpolated to DARDAR grid.
       
        Returns
        -------
        grid_t2m : np.array containing the interpolated values
        dimensions [lat, lon]

        """
        var           = "2m_temperature"
        ERA_t2        = ERA5s(self.t_0, self.t_1, var)
        grid_t2m       = ERA_t2.interpolate(self.dardar)
        grid_t2m      = np.expand_dims(grid_t2m, axis = 1)
        
        return grid_t2m

    @property
    def wind_speed(self):
        """
        ERA5 10m u and v fields interpolated to DARDAR grid.
        The wind speed is calculated from u and v vectors
       
        Returns
        -------
        wind_speed : np.array containing the interpolated values
        dimensions [lat, lon]

        """
        # u 10m 
        var           = "10m_u_component_of_wind"
        ERA_u         = ERA5s(self.t_0, self.t_1, var)
        grid_u        = ERA_u.interpolate(self.dardar)
        
        # v 10 m 
        var           = "10m_v_component_of_wind"
        ERA_v         = ERA5s(self.t_0, self.t_1, var)
        grid_v        = ERA_v.interpolate(self.dardar)
        
        # convert to wind speed
        wind_speed    = np.sqrt(grid_u**2 + grid_v**2)
        wind_speed    = np.expand_dims(wind_speed, axis = 1)
        
        return wind_speed
    
    @property    
    def wind_direction(self):
        """
        ERA5 10m u and v fields interpolated to DARDAR grid.
        The wind speed is calculated from u and v vectors
       
        Returns
        -------
        wind_direction : np.array containing the interpolated values
        dimensions [lat, lon]

        """
        
        # u 10m 
        var           = "10m_u_component_of_wind"
        ERA_u         = ERA5s(self.t_0, self.t_1, var)
        grid_u        = ERA_u.interpolate(self.dardar)
        
        # v 10 m 
        var           = "10m_v_component_of_wind"
        ERA_v         = ERA5s(self.t_0, self.t_1, var)
        grid_v        = ERA_v.interpolate(self.dardar)
        
        # convert to wind direction
        wind_dir      = np.arctan2(grid_u, grid_v) * 180 / np.pi  # degree
        
        wind_dir    = np.expand_dims(wind_dir, axis = 1)
        
        return wind_dir
    
    @property
    def vmr_h2o(self, p_grid = None):   
        """
        interpolated ERA5 VMR fields to DARDAR grid and pressure grid
        defined in self.p_grid
        
        VMR is calculated using specific humidity. The conversion is simply
        
        vmr = Q/1-Q
        
        The VMR is firstly calculated on ERA5 pressure grid and later 
        interpolated to desired p_grid as log(vmr)
        
        The code for conversion of relative humidity is commented. 
        The RH2VMR conversion assumes the defintions of hydrostatic
        equilibrium from IFS

        Returns
        -------
        grid_q2vmr : np.array containing the interpolated values
        dimensions [1, p, lat, lon]

        """
        
        p           = self.p_grid   
        
        # var             = "relative_humidity"
        # ERA_r           = ERA5p(t_0, t_1, var)  # loads ERA5 data as xarray 
        # grid_r          = ERA_r.interpolate(dardar, p_grid = None) 
        
        # if np.any(grid_r < 0):
        #     print ("encountered negative values")
         
        # temperature
            
        # p_grid          = ERA_q.era["level"].data * 100 # Pa
        # grid_p          = np.tile(p_grid, (grid_t.shape[1], 1))
        # grid_p          = grid_p.T 
        # r2vmr           = rh2vmr(grid_p , grid_t, grid_r)
        
        
#        grid_t          = self.temperature(p_grid = None)
            
#       interpolate log(vmr) as a function of log(p)
        
        # f_r             =  interpolate.interp1d(np.log(p_grid), np.log(r2vmr), axis = 0 )
        # grid_r2vmr      =  np.exp(f_r(np.log(p)))
        
        var             = "specific_humidity"
        ERA_q           = ERA5p(self.t_0, self.t_1, var)
        grid_q          = ERA_q.interpolate(self.dardar, p_grid = None)
        q2vmr           = thermodynamics.specific_humidity2vmr(grid_q)    
        p_era           = ERA_q.era["level"].data * 100 # Pa
            
        f_q             =  interpolate.interp1d(np.log(p_era), np.log(q2vmr), axis = 0 )
        grid_q2vmr      =  np.exp(f_q(np.log(p* 100)))
        
        grid_q2vmr      = np.expand_dims(grid_q2vmr, axis = (0, 3))
        
        return grid_q2vmr
    
    @property
    def vmr_N2(self):
        """
        VMR values for N2.
        The standard value for N2 is scaled by VMR of h2o
        Returns
        -------
        grid_N2 : np.array containing the interpolated values in
        dimensions [1, p, lat, lon]
        """        
        grid_q2vmr       = self.vmr_h2o
        
        grid_N2         = np.ones(grid_q2vmr.shape) * 0.781 # vmr N2
        
        grid_N2         = scale_vmr(grid_N2, grid_q2vmr)
        
        return grid_N2

        
    @property    
    def vmr_O2(self): 
        """
        VMR values for O2.
        The standard value for O2 is scaled by VMR of h2o
        Returns
        -------
        grid_N2 : np.array containing the interpolated values in
        dimensions [1, p, lat, lon]
        """ 
        
        grid_q2vmr       = self.vmr_h2o
        
        
        grid_O2         = np.ones(grid_q2vmr.shape) * 0.209 # vmr O2

        grid_O2         = scale_vmr(grid_O2, grid_q2vmr)
        
        return grid_O2
        
    @property
    def iwc(self):
        """
        The IWC data from DARDAR interpolated to pressure grid defined in
        p_grid
        -------
        grid_N2 : np.array containing the interpolated values in
        dimensions [1, p, lat, lon]
        """ 
        
        p               = self.p_grid   
        
        iwc             = self.dardar.get_data('iwc')
        height_d        = self.dardar.get_data('height')
        p_grid_d        = alt2pres(height_d)
        f               = interpolate.interp1d(p_grid_d, iwc)
        grid_iwc        = f(p * 100)
        grid_iwc        = grid_iwc.T
        grid_iwc        = np.expand_dims(grid_iwc, axis = (0, 3))
        grid_iwc        = np.expand_dims(grid_iwc, 3)   
        
        return grid_iwc                    
        
            

        
        
        