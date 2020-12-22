#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 09:43:17 2020

@author: inderpreet
"""
import os
import numpy as np
from datetime import datetime, timedelta
from era2dardar.DARDAR import DARDARProduct
from era2dardar.ERA5 import ERA5p, ERA5s
from era2dardar.utils.alt2pressure import alt2pres, pres2alt
import typhon.physics.thermodynamics as thermodynamics 
import typhon.physics.atmosphere as atmosphere
from era2dardar.utils.scale_vmr import scale_vmr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from era2dardar.utils.thermodynamics  import mixr2massconc, airdensity, vmr2mixr, rh2vmr
from era2dardar.utils.thermodynamics import massconc2tcwv
from integrate_vmr import integrate_vmr
from scipy import interpolate
import typhon.arts.xml as xml
from filter_LWC import filter_LWC
from Z2dbZ import Z2dbZ
plt.rcParams.update({'font.size': 18})

#%%--------------------------------------------------------------------------------------
# DARDAR file
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310114257_50670.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071190249_47194.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071104824_47189.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/06/07/DARDAR-CLOUD_v2.1.1_2015158155052_48459.hdf")

filename = file

# create DARDAR instance
dardar = DARDARProduct(filename, latlims = [0, 10])
dardar.plot_overpass()  
# time stamp of DARDAR data 
t_0 = dardar.filename2date()
t_1 = t_0 + timedelta(minutes = 30)

#%%
# final pressure grid 
z_grid          = np.arange(-100, 25000, 250) 
p               = alt2pres(z_grid) # final pressure grid [Pa]


#%%
#    read ERA5 Q, RH and T on ERA5 grid, and later interpolate to final grid
#    as vmr(log(p))  
    
   
var             = "relative_humidity"
ERA_r           = ERA5p(t_0, t_1, var)  # loads ERA5 data as xarray 
grid_r          = ERA_r.interpolate(dardar, p_grid = None) 

if np.any(grid_r < 0):
    print ("encountered negative values")
 
# temperature
var             = "temperature"
ERA_t           = ERA5p(t_0, t_1, var)
grid_t          = ERA_t.interpolate(dardar, p_grid = None)

var             = "specific_humidity"
ERA_q           = ERA5p(t_0, t_1, var)
grid_q          = ERA_q.interpolate(dardar, p_grid = None)
q2vmr           = thermodynamics.specific_humidity2vmr(grid_q)


p_grid          = ERA_q.era["level"].data * 100 # Pa
grid_p          = np.tile(p_grid, (grid_t.shape[1], 1))
grid_p          = grid_p.T 
r2vmr           = rh2vmr(grid_p , grid_t, grid_r)


#interpolate log(vmr) as a function of log(p)

f_r             =  interpolate.interp1d(np.log(p_grid), np.log(r2vmr), axis = 0 )
grid_r2vmr      =  np.exp(f_r(np.log(p)))

f_q             =  interpolate.interp1d(np.log(p_grid), np.log(q2vmr), axis = 0 )
grid_q2vmr      =  np.exp(f_q(np.log(p)))

#%%
# calculate vmr for N2 and O2
# the standard values are scaled by vmr of h20

grid_N2         = np.ones(grid_q2vmr.shape) * 0.781 # vmr N2
grid_O2         = np.ones(grid_q2vmr.shape) * 0.209 # vmr O2

grid_N2         = scale_vmr(grid_N2, grid_q2vmr)
grid_O2         = scale_vmr(grid_O2, grid_q2vmr)

#%%
#  interpolate LWC on final pressure grid
var             = "temperature"
ERA_t           = ERA5p(t_0, t_1, var)
grid_t          = ERA_t.interpolate(dardar, p_grid = p * 0.01)

# volume mixing ratio
var             = "specific_cloud_liquid_water_content"
ERA_lwc         = ERA5p(t_0, t_1, var)
grid_lwc        = ERA_lwc.interpolate(dardar, p_grid = p * 0.01)

grid_p          = np.tile(p, (grid_lwc.shape[1], 1))

rho             = thermodynamics.density(grid_p.T, grid_t) # density of air
grid_lwc        = grid_lwc * rho  # convert lwc units to mass concentration    

# Z
Z               = dardar.get_data('Z')
height_d        = dardar.get_data('height')
p_grid_d        = alt2pres(height_d)
f               = interpolate.interp1d(p_grid_d, Z)
grid_Z          = f(p)
grid_Z          = grid_Z.T
mask            = grid_Z == 0

dbZ             = grid_Z.copy()
dbZ[~mask]      = Z2dbZ(grid_Z[~mask])
grid_lwc        = filter_LWC(grid_lwc, dbZ)

mask_dbZ         = dbZ <= -25.0
grid_Z[mask_dbZ] = 0
grid_Z           = np.expand_dims(grid_Z, 2)

# IWC
iwc             = dardar.get_data('iwc')
height_d        = dardar.get_data('height')
p_grid_d        = alt2pres(height_d)
f               = interpolate.interp1d(p_grid_d, iwc)
grid_iwc        = f(p)
grid_iwc        = grid_iwc.T
grid_iwc        = np.expand_dims(grid_iwc, 0)
grid_iwc        = np.expand_dims(grid_iwc, 3)

#%%
vmr_grid    = np.stack([grid_N2,
                      grid_O2,
                      grid_q2vmr,
                      grid_lwc])



xml.save(vmr_grid, 'xmls/vmr_field.xml')
xml.save(grid_Z, 'xmls/reflectivities.xml')
xml.save(grid_iwc, 'xmls/IWC.xml')

#%%
# plot the interpolated q2vmr and r2vmr fields
 
p_grid = p * 0.01  ## for plotting use hPa
lat_d    = dardar.get_data('latitude')

# contour plots for VMR
fig, axs = plt.subplots(3, 1, figsize = [20, 25])

fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2, hspace = 0.5)

#for var, ax in zip([r2vmr, q2vmr], axs.ravel()):
for var, ax in zip([grid_r2vmr, grid_q2vmr], axs.ravel()):
    im = (ax.pcolormesh(lat_d, p_grid, var, shading = 'auto', cmap = 'tab20c',
#                        vmin = 0.000001, vmax = 0.025, 
                        norm=colors.LogNorm() ))
    
    fig.colorbar(im, ax=ax, extend = 'both')

    ax.set_xlabel('Latitude [deg]')
    ax.set_ylabel('Pressure [hPa]')

    ax.set_ylim(ax.get_ylim()[::-1])

    ax.set_yscale('log')

axs[0].set_title('RH to VMR')
axs[1].set_title('Q to VMR')
axs[2].set_title('(r2vmr - q2vmr)/r2vmr')

d = ((grid_q2vmr - grid_r2vmr)/grid_r2vmr * 100)
#d[t_mix] = np.nan
#d[t_liq] = np.nan
im = (axs[2].pcolormesh(lat_d, p_grid, d, shading = 'auto', cmap = 'coolwarm',
                        vmin = -20, vmax = 20))

fig.colorbar(im, ax=axs[2], extend = 'both')
axs[2].set_xlabel('Latitude [deg]')
axs[2].set_ylabel('Pressure [hPa]')
axs[2].set_yscale('log')
axs[2].set_ylim(axs[2].get_ylim()[::-1])

fig.savefig('Figures/interpolated_vmr.png', bbox_inches = 'tight')





