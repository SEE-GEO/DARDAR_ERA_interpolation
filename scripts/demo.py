#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:00:30 2020

@author: inderpreet
"""
import os
import glob
import numpy as np
from datetime import datetime, timedelta
from era2dardar.DARDAR import DARDARProduct
from era2dardar.ERA5 import ERA5p, ERA5s
from era2dardar.atmData import atmdata
from era2dardar.utils.alt2pressure import alt2pres
import typhon.physics.thermodynamics as thermodynamics 
import typhon.physics.atmosphere as atmosphere
import typhon.arts.xml as xml
from typhon.topography import SRTM30
from era2dardar.utils.scale_vmr import scale_vmr
from era2dardar.dardar2atmdata import dardar2atmdata
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import shutil


# pressure grid
p_grid = alt2pres(np.arange(-100, 25000, 250))

# latitudinal extent
latlims    = [-30, 30]

# year and month of data
year = "2009"
month = "07"

# add all eligible files to dardarfiles
inpath = os.path.join(os.path.expanduser("~/Dendrite/SatData/DARDAR"), year, month)
dardarfiles = glob.glob(os.path.join(inpath, "*", "*.hdf"))

outpath = os.path.expanduser("~/Dendrite/Projects/IWP/")


# DARDAR file
# file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310114257_50670.hdf")
# #file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071190249_47194.hdf")
# #file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071104824_47189.hdf")
# #file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/06/07/DARDAR-CLOUD_v2.1.1_2015158155052_48459.hdf")

# file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/06/03/DARDAR-CLOUD_v2.1.1_2015154193320_48403.hdf")

file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2009/06/03/DARDAR-CLOUD_v2.1.1_2009154021055_16481.hdf")


# start the loop

for dardarfile in dardarfiles[10:20]:
    for N in ["A","D"]:
 
        try:
            dardar = DARDARProduct(dardarfile, latlims = latlims, node = N)
            dardar.plot_scene()
        except:
            print ("descending pass not available")
            continue
            
# get all atmfields
        atm_fields  = dardar2atmdata(dardar, p_grid)
        
        date = dardar.filename2date()
        outdir = year + "_" + date.strftime("%3j") + "_" + date.strftime("%2H") + "_" + N      
        
        outdir = os.path.join(outpath, outdir)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        
  
# save xml files          
        for key in atm_fields.keys():
            
             filename = os.path.join(outdir,   key + ".xml")
            
             xml.save(atm_fields[key], filename)

        output_filename = os.path.basename(outdir) 
        shutil.make_archive(outdir, 'zip', outdir)      
 
# remove downloaded ERA files to avoid memory build up     
    erafiles = glob.glob("ERA5/*/*")
    for f in erafiles:    
        os.remove(f)        
        
    
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

