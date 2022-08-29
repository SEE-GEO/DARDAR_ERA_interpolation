#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:23:46 2020

Dictionary containing long- and short-names of variables available 
in ERA5 pressure or surface data.

More information on the 
single level variables can be found :
https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview

and 

pressure levels can be found :
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview
    

@author: inderpreet
"""



parameters = {"temperature": "t",
                     "geopotential": "z",
                     "specific_humidity" : "q",
                     "relative_humidity" : "r",
                     "specific_cloud_ice_water_content": "ciwc",
                     "specific_cloud_liquid_water_content": "clwc",
                     "10m_v_component_of_wind": "v10",
                     "10m_u_component_of_wind": "u10",
                     "v_component_of_wind": "v",
                     "u_component_of_wind": "u",
                     "skin_temperature": "skt",
                     "2m_temperature": "t2m",
                     "total_column_water_vapour": 'tcwv',
                     "surface_pressure": "sp",
                     "sea_ice_cover": "siconc",
                     "land_sea_mask": "lsm",
                     "snow_depth": "sd",
                     "ozone_mass_mixing_ratio": "o3",
                     "surface_pressure": "sp",
                     "orography": "z"}
