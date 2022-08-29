#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:27:55 2021

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
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import shutil
import matplotlib.pyplot as plt

from era2dardar.atmData import atmdata

# pressure grid
p_grid = alt2pres(np.arange(-700, 20000, 250))
p_grid = (np.concatenate([p_grid, 
                         np.array([30, 20, 10, 7, 5, 3, 2, 1]) * 100]))

#dardarfile = "/home/inderpreet/Dendrite/SatData/DARDAR/2010/01/31/DARDAR-CLOUD_v2.1.1_2010031065234_20008.hdf" 

dardarfile = "/home/inderpreet/Dendrite/SatData/DARDAR/2010/01/31/DARDAR-CLOUD_v2.1.1_2010031033448_20006.hdf"
#dardarfile = "/home/inderpreet/Dendrite/SatData/DARDAR/2010/01/20/DARDAR-CLOUD_v2.1.1_2010020035329_19846.hdf"

#dardarfile = "/home/inderpreet/Dendrite/SatData/DARDAR/2009/07/31/DARDAR-CLOUD_v2.1.1_2009212191541_17336.hdf"
#dardarfile = "/home/inderpreet/Dendrite/SatData/DARDAR/2010/01/29/DARDAR-CLOUD_v2.1.1_2010029201601_19987.hdf"
#dardarfile = "/home/inderpreet/Dendrite/SatData/DARDAR/2009/12/01/DARDAR-CLOUD_v2.1.1_2009335203452_19128.hdf"

# z_surface = xml.load("/home/inderpreet/data/temp/z_surface.xml")

# z_field = xml.load("/home/inderpreet/data/temp/z_field.xml")

latlims = [-65, 65]
N = "A"
dardar = DARDARProduct(dardarfile, latlims = latlims, node = N)
dardar.plot_scene()

print ("t_0, t_1", dardar.t_0, dardar.t_1)


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



atm = atmdata(dardar, p_grid, domain)

z_surface_new = atm.z_surface
z_field_new      = atm.z_field
lat = atm.lat
z_old = atm.z_field_old


fig, ax  = plt.subplots(1, 1, figsize = [8, 8])

ax.plot(lat, z_field_new[21, :, 0], label = "new")
ax.plot(lat, z_old[21, :, 0], label = "old")
ax.plot(lat, z_surface_new, label = "surface")
ax.set_ylabel("altitude[m]")
ax.set_xlabel("latitude")
ax.legend()
fig.savefig("ERA5_elevation.png", bbox_inches = "tight")

fig, ax  = plt.subplots(1, 1, figsize = [8, 8])
#ax2  =ax.twiny()
diff = z_old - z_field_new
ax.plot(lat, diff[21, :, 0], label = "z_field")
#ax2.plot(lat, z_surface_new[:], color = "red", label = "surface")
ax.set_ylabel("difference altitude 500 hPa [m]")
ax.set_xlabel("latitude")
ax.legend()
#ax2.legend()
fig.savefig("diff_z_field.png", bbox_inches = "tight")


for i in range(len(p_grid)):
    print (np.abs(np.diff(z_field_new[i, :, 0])).max())

 
