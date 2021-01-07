#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 26 12:33:27 2020

@author: inderpreet
"""

from pyhdf.SD import SD, SDC
import os
from datetime import datetime, timedelta
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from era2dardar.utils.seconds2datetime import seconds2datetime
from era2dardar.utils.Z2dbZ import Z2dbZ


class DARDARProduct():
    """
    The DARDARProduct class defines an interface to read DARDAR profiles
    
    """    
    
    def __init__(self, filename, latlims = None, node = "A"):
        """
        Opens the dardar hdf4 dataset

        Parameters
        ----------
        filename  (str): input DARDAR file
        latlims (None, list)= None, global data is returned; 
                [lat1, lat2], list containing lower and upper limits of
        latitude values used to subset DARDAR pass
        
        node : string "A" or "D_N" or "D_S", the ascending node or descending node(NH or SH)
               The default node selected is ascending one
               If Latlims is None, the complete scene is used.
               
        Exceptions:
            
            No data is returned if the latitude limits are out of scene


        """
        self.filename = filename
#        file_name = 'MYD06_L2.A2007219.2010.006.2014053202546.hdf'
        self.file = SD(self.filename, SDC.READ)
        
        datasets_dic = self.file.datasets()

# list of SDS variables
        self.SDS = datasets_dic.keys()
        
        self.node = node
             
        self.latlims = latlims
        
        if self.latitude.size == 0:
            raise Exception("No data returned, input another latlims")
    
    @property
    def latitude(self):
        """
        gets latitude values for DARDAR pass

        Returns
        -------
        lat : np.array containing latitudes [deg]

        """
        
        lat = self.get_data("latitude")
        return lat
    
    @property
    def longitude(self):
        """
        gets longitude values for DARDAR pass

        Returns
        -------
        lon : np.array containing longitudes [deg]

        """     
        lon = self.get_data("longitude")
        return lon

    @property  
    def iwc(self):
        """
        gets IWC values for DARDAR pass

        Returns
        -------
        iwc : np.array containing ice water content []

        """
        
        iwc = self.get_data("iwc")
        return iwc
    
    @property    
    def Z(self):
        """
        gets Z values for DARDAR pass

        Returns
        -------
        z : np.array containing reflectivities []

        """
        Z = self.get_data("Z")
        return Z
    
    @property
    def time(self):
        """
        gets time values for DARDAR pass

        Returns
        -------
        z : np.array containing datetime objects for time

        """
        time = self.get_data("time")
        return time
    
    @property
    def t_0(self):
        """
        get the datetime object for first DARDAR profile in the selected scene
        converts seconds to datetime object using method seconds2datetime
        
        Returns
        -------
        datetime object
        """
        
        t = self.time[0]
        date = self.filename2date()
        date = date.replace(hour = 0, minute = 0)
        t = date + timedelta(seconds = float(t))
        return t
    
    @property
    def t_1(self):
        """
        get the datetime object for last DARDAR profile in the selected scene
        converts seconds to datetime object using method seconds2datetime
        
        Returns
        -------
        datetime object
        """
        
        t = self.t_0 + timedelta(minutes = 30)
        return t
        
    
    @property
    def height(self):
        """
        gets altitudes over which DARDAR profiles are defined

        Returns
        -------
        height : np.array containing altitudes [m]

        """
        height = self.get_data("height")
        return height
        
        
    def get_data(self, variable):
        """
        get the data for the selected variable
        complete list of SDS variables is also contained in self.SDS

        Parameters
        ----------
        variable : Input SDS variable,

        Returns
        -------
        ndarray containing the input SDS variable

        """
        if variable not in self.SDS:
            raise Exception("'Valid SDS should be one of ", self.SDS) 
            
        sds_obj = self.file.select(variable) # select sds

        data = sds_obj.get() # get sds data    
            
        if self.latlims is not None:
            # subsetting  data   
            if variable != "height":
                lat1, lat2 = self.latlims
                lat  = self.file.select('latitude').get()
                dn_flag  = self.file.select("day_night_flag").get()
            
                inds = np.where((lat >= lat1) & (lat <= lat2))[0]
                
                if len(inds) == 0:
                    raise Exception("No data in the latitude limits given, please check limits again")
                lat_sub = lat[inds]

# to avoid extracting two latitudes, each from ascending and descending node  
            
                diff = (np.diff(lat_sub,
                        append = lat_sub[-1] + lat_sub[-1] - lat_sub[-2]))

                        
                mask1 = diff > 0               
                
                if np.all(~mask1):
                    raise Exception("No increasing latitudes were found")
                    
# day night mask
                mask2 = dn_flag[inds] == 0     

                if self.node == "A": 
                    mask = np.logical_and(mask1, mask2)                    
                    data = data[inds][mask]
                if self.node == "D_S":
                    mask3 = lat_sub < 0
                    mask  = np.logical_and(~mask1, ~mask2) 
                    mask  = np.logical_and(mask, mask3)
                    if np.sum(mask) == 0:
                        raise Exception("No data in the SH descending pass")
                    data  = data[inds][mask]        
                if self.node == "D_N":
                    mask3 = lat_sub >= 0
                    mask  = np.logical_and(~mask1, ~mask2)  
                    mask  = np.logical_and(mask, mask3)
                    if np.sum(mask) == 0:
                        raise Exception("No data in the NH descending pass")
                    data  = data[inds][mask]                      
                
                
        return data
    
   
    def filename2date(self):
        """
        extracts time stamp from dardar filename

        Returns
        -------
        datetime object with the time stamp of DARDAR data
        """
        
        filename = os.path.basename(self.filename)
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
    
    
    def plot_scene(self):
        """
        plots the overpass of DARDAR

        Returns
        -------
        None.

        """
        plt.figure(figsize=(12, 6))
        m = Basemap(llcrnrlon=0.,llcrnrlat=-85.,urcrnrlon=360.,urcrnrlat=85.,\
                      rsphere=(6378137.00,6356752.3142),\
                      resolution='c',projection='cyl')
        m.shadedrelief(scale = 0.1)
        
        lon_d = self.get_data("longitude")
        lat_d = self.get_data("latitude")
        lon_d = lon_d % 360
        m.scatter(lon_d[:],lat_d[:], latlon = True)
        plt.savefig('dardar_pass.png', bbox_inches = 'tight')
        plt.show()
        