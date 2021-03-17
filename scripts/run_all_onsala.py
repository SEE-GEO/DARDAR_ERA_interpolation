#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:00:30 2020

interface to generate atmospheric data for ONSALA

@author: inderpreet
"""
import os
import glob
import numpy as np
import random

from era2dardar.onsala_atmdata import onsala_atmdata 
from era2dardar.locations import locations
from era2dardar.utils.alt2pressure import alt2pres
import typhon.arts.xml as xml
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import shutil



def run_all_cases(p_grid, onsala, t0, outpath): 

  
    
        print ("doing time ", t0)

        outdir = t0.strftime("%Y") + "_" + t0.strftime("%3j") + "_" + t0.strftime("%2H")    
        
        print (outdir)
        if os.path.isfile(os.path.join(outpath, outdir + '.zip')):
            print ('file %s already exists, doing next file'%outdir)
            t0 +=  timedelta(days = 1)
            return
                
        
        print ("t_0, t_1", onsala.t_0, onsala.t_1)

        
        lon          = onsala.longitude 
        lat          = onsala.latitude
        
        lat1         = np.around(lat.min() - 1)
        lat2         = np.around(lat.max() + 1)
        
        lon1         = np.around(lon.min() - 1)
        lon2         = np.around(lon.max() + 1)

        if lon1 < -180:
            lon1 = -180.0
        if lon2 > 180:
            lon2 = 180.0    
  
        # when encountering prime meridian, download global data
        # keeps interpolation simple          
        if lon1 == -180 or lon2 == 180:
            lon1 = -180.
            lon2 =  180.
            
        # domain for which ERA5 data is downloaded           
        domain  = [lat1, lat2, lon1, lon2]
        
        # get all atmfields as a directory
        atm_fields  = onsala_atmdata(onsala, p_grid, domain = domain)
        
        
        outdir = os.path.join(outpath, outdir)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        
  
        # save xml files  to a zipped folder        
        for key in atm_fields.keys():
            
             filename = os.path.join(outdir,   key + ".xml")
            
             xml.save(atm_fields[key], filename)


        
        shutil.make_archive(outdir, 'zip', outdir)   

        # remove unzipped folder
        shutil.rmtree(outdir)
        
 
        # remove downloaded ERA files  
        erafiles = (glob.glob(os.path.join("ERA5/*/", "*" 
                                           + t0.strftime("%4Y") 
                                           + t0.strftime("%2m") 
                                           + t0.strftime("%2d")  
                                           + t0.strftime("%2H") +"*")))
        for f in erafiles:    
            os.remove(f)    

        

if __name__ == "__main__":

    # pressure grid
    p_grid = alt2pres(np.arange(-700, 20000, 250))
    p_grid = (np.concatenate([p_grid, 
                             np.array([30, 20, 10, 7, 5, 3, 2, 1]) * 100]))

    
    t0 = datetime(2013, 1, 1, 12, 00, 00 )
    t1 = datetime(2014, 1, 1, 12, 00, 00)
    
    # Onsala location
    lat = 57. + 23./60. + 35./3600.
    lon = 11. + 55./60. + 04./3600.
        
    while t0 < t1:
        
        onsala = locations(lat, lon, t0)
    
        outpath = os.path.expanduser("~/Dendrite/Projects/Onsala_ERA_data/")    
    
        run_all_cases(p_grid, onsala, t0, outpath)
        
        t0 +=  timedelta(days = 1)    