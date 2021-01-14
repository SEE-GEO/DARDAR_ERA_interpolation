#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:00:30 2020

interface to generate xml files for multiple DARDARfiles

@author: inderpreet
"""
import os
import glob
import numpy as np
import random

from era2dardar.dardar2atmdata import dardar2atmdata
from era2dardar.DARDAR import DARDARProduct
from era2dardar.utils.alt2pressure import alt2pres
import typhon.arts.xml as xml

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import shutil



def run_all_cases(p_grid, dardarfiles, Nodes, latlims, inpath, outpath, year, month):
    
    for dardarfile in dardarfiles:
        
        for N in Nodes:
     
            try:
                dardar = DARDARProduct(dardarfile, latlims = latlims, node = N)
                dardar.plot_scene()
            except:
                raise Exception("descending pass not available")
            
            print ("t_0, t_1", dardar.t_0, dardar.t_1)
    
            date = dardar.filename2date()
            outdir = year + "_" + date.strftime("%3j") + "_" + date.strftime("%2H") + "_" + N      
            
    
# check if file already exists
            if os.path.isfile(os.path.join(outpath, outdir + '.zip')):
                print ('file %s already exists, doing next file'%outdir)
                continue
            
    
            
            lon          = dardar.longitude 
            lat          = dardar.latitude
            
            lat1         = np.around(lat.min() - 2)
            lat2         = np.around(lat.max() + 2)
            
            lon1         = np.around(lon.min() - 2)
            lon2         = np.around(lon.max() + 2)
    
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
            atm_fields  = dardar2atmdata(dardar, p_grid, domain = domain)
            
#        date = dardar.filename2date()
#        outdir = year + "_" + date.strftime("%3j") + "_" + date.strftime("%2H") + "_" + N      
            
            outdir = os.path.join(outpath, outdir)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            
      
# save xml files  to a zipped folder        
            for key in atm_fields.keys():
                
                 filename = os.path.join(outdir,   key + ".xml")
                
                 xml.save(atm_fields[key], filename)
    
            output_filename = os.path.basename(outdir) 
            shutil.make_archive(outdir, 'zip', outdir)   
#        shutil.make_archive(outdir, 'zip', 'arts.xml') 
            shutil.rmtree(outdir)
            
     
# remove downloaded ERA files  
        erafiles = glob.glob("ERA5/*/*")
        for f in erafiles:    
            os.remove(f)        
        

if __name__ == "__main__":

    # pressure grid
    p_grid = alt2pres(np.arange(-700, 20000, 250))
    p_grid = np.concatenate([p_grid, np.array([30, 20, 10, 7, 5, 3, 2, 1]) * 100])
    
    # latitudinal extent
    latlims    = [-30, 30]
    
    # year and month of data
    year = "2009"
    month = "08"
    
    # add all eligible files to dardarfiles
    inpath = os.path.join(os.path.expanduser("~/Dendrite/SatData/DARDAR"), year, month)
    dardarfiles = glob.glob(os.path.join(inpath, "*", "*.hdf"))
    
    outpath = os.path.expanduser("~/Dendrite/Projects/IWP/IceCube/DARDAR_ERA_m30_p30")
    
    Nodes =  ["A", "D_S", "D_N"]
    # random shuffle dardarfiles
    #random.shuffle(dardarfiles)
    
    # start the loop for all cases
    
    run_all_cases(p_grid, dardarfiles, Nodes, latlims, inpath, outpath, year, month)
    
    
        
    #------------------------------------to check data----------------------------    
    # dardar.plot_scene()
    
    # #check IWC values
    # fig, ax = plt.subplots(1,1, figsize = [16, 8])
    # cmap = 'coolwarm'
    # lat_d = dardar.latitude
    # height_d  = dardar.height
    # iwc       = dardar.iwc
    # p         = alt2pres(height_d) * 0.01
    # im = ax.pcolormesh(lat_d, p, iwc.T, cmap=cmap, norm=colors.LogNorm(),
    #                                       vmin = 0.000001, vmax = 0.01)
    # ax.set_ylim(ax.get_ylim()[::-1])
    # # #ax.set_yscale('log')
    # fig.colorbar(im, ax=ax, label = 't [K]' , extend = 'both')
    # #
    # # get atmdata on DARDAR grid
    
    # atm = atmdata(dardar, p_grid)
    
    # #lsm = atm.sea_ice_cover
    # #plt.plot(lat_d, lsm)
    # z = atm_fields["z_field"]
    # fig, ax = plt.subplots(1,1, figsize = [8, 8])
    # ax.plot(np.log(p_grid * 0.01), z[:, 1000, 0]/1000)
    # ax.plot(np.log(p_grid * 0.01), z[:, 2000, 0]/1000)
    
    # ax.plot(np.log(p_grid * 0.01), z[:, 3, 0]/1000)
