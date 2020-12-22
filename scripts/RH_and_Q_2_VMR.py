#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 12:33:03 2020

comparing the VMR conversion from RH and Q

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
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from era2dardar.utils.thermodynamics  import mixr2massconc, airdensity, vmr2mixr, rh2vmr
from era2dardar.utils.thermodynamics import massconc2tcwv
from integrate_vmr import integrate_vmr
from scipy import interpolate
plt.rcParams.update({'font.size': 18})


#%%--------------------------------------------------------------------------------------
# DARDAR file
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310114257_50670.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071190249_47194.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071104824_47189.hdf")
#file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/06/07/DARDAR-CLOUD_v2.1.1_2015158155052_48459.hdf")

filename = file

# create DARDAR instance
dardar = DARDARProduct(filename,latlims = None)
dardar.plot_overpass()  
# time stamp of DARDAR data 
t_0 = dardar.filename2date()
t_1 = t_0 + timedelta(minutes = 30)

#%%-----------------------------------------------------------------------------------------
# interpolate required ERA data
p_grid = None

var             = "temperature"
ERA_t           = ERA5p(t_0, t_1, var)
grid_t          = ERA_t.interpolate(dardar, p_grid)

# pressure grid
z_grid          = np.arange(-100, 25000, 250) 
p_grid          = ERA_t.era["level"].data[3:] * 100 # [Pa]
#p_grid         = alt2pres(z_grid) 


var             = "surface_pressure"
ERA_sp          = ERA5s(t_0, t_1, var)
grid_sp         = ERA_sp.interpolate(dardar)

var             = "total_column_water_vapour"
ERA_tcwv        = ERA5s(t_0, t_1, var)
grid_tcwv       = ERA_tcwv.interpolate(dardar)

var             = "specific_humidity"
ERA_q           = ERA5p(t_0, t_1, var)
grid_q          = ERA_q.interpolate(dardar, p_grid * 0.01)
q2vmr           = thermodynamics.specific_humidity2vmr(grid_q)



var             = "relative_humidity"
ERA_r           = ERA5p(t_0, t_1, var)
grid_r          = ERA_r.interpolate(dardar, p_grid * 0.01) 

if np.any(grid_r < 0):
    print ("encountered negative values")

# temperature
var             = "temperature"
ERA_t           = ERA5p(t_0, t_1, var)
grid_t          = ERA_t.interpolate(dardar, p_grid * 0.01)


#%% ------------------convert relative humidity to vmr --------------------------
grid_p          = np.tile(p_grid, (grid_t.shape[1], 1))
grid_p          = grid_p.T 
r2vmr           = rh2vmr(grid_p , grid_t, grid_r)


#%%--------------------------------------------------------------------------------
#interpolate log(vmr) as a function of log(p)

# f_r       =  interpolate.interp1d(np.log(p_grid), np.log(r2vmr), axis = 0 )
# r2vmr     =  np.exp(f_r(np.log(p_grid)))

# f_q       =  interpolate.interp1d(np.log(p_grid), np.log(q2vmr), axis = 0 )
# q2vmr     =  np.exp(f_q(np.log(p_grid)))


#%% -----------------------------------------------------------------------------------
# convert vmr to mass concentration
M_w         = 18.0153e-3
mixr_r      = vmr2mixr(r2vmr, M_w)
mixr_q      = vmr2mixr(q2vmr, M_w)

mc_r        = mixr2massconc(mixr_r, grid_p , grid_t)
mc_q        = mixr2massconc(mixr_q, grid_p , grid_t)

#%% -----------------------------------------------------------------------------------
##  integrate mass concentration to get TCWV
  
tcwv_r   = massconc2tcwv(mc_r, p_grid, grid_sp)    
tcwv_q   = massconc2tcwv(mc_q, p_grid, grid_sp)   

fig, ax = plt.subplots(1, 1, figsize= [8, 8])

lat_d = dardar.get_data('latitude')    

ax.plot(grid_tcwv, lat_d, label = 'ERA5 TCWV ')
ax.plot( tcwv_q, lat_d, label = 'Q to TCWV')
ax.plot( tcwv_r, lat_d,label = 'RH to TCWV')
ax.set_ylabel('Latitude [deg]')
ax.set_xlabel('TCWV [kg/m2]')    
ax.legend()
fig.savefig('Figures/comparison_TWCV.png', bbox_inches = 'tight')

    
#%%-------------------# plotting data below----------------------------------------

p_grid = p_grid * 0.01  ## for plotting use hPa
lat_d    = dardar.get_data('latitude')

# contour plots for VMR
fig, axs = plt.subplots(3, 1, figsize = [20, 25])

fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2, hspace = 0.5)

#for var, ax in zip([r2vmr, q2vmr], axs.ravel()):
for var, ax in zip([r2vmr, q2vmr], axs.ravel()):
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

d = ((q2vmr - r2vmr)/r2vmr * 100)
#d = ((r2vmr - np.mean(r2vmr, axis = 1).reshape(-1, 1)))
#d[t_mix] = np.nan
#d[t_liq] = np.nan
im = (axs[2].pcolormesh(lat_d, p_grid, d, shading = 'auto', cmap = 'coolwarm',
                        vmin = -20, vmax = 20,
                        ))

fig.colorbar(im, ax=axs[2], extend = 'both')
axs[2].set_xlabel('Latitude [deg]')
axs[2].set_ylabel('Pressure [hPa]')
axs[2].set_yscale('log')
axs[2].set_ylim(axs[2].get_ylim()[::-1])

# t = grid_t.copy()
# t[t_mix] = np.nan


# t_mean = np.mean(grid_t, axis = 1).reshape(-1, 1)
# im = (axs[3].pcolormesh(lat_d, p_grid, grid_t, shading = 'auto', cmap = 'coolwarm', ))
# axs[3].set_ylim(axs[3].get_ylim()[::-1])
# fig.colorbar(im, ax=axs[3], extend = 'both')
# axs[3].set_yscale('log')
# axs[3].set_title('temperature')
# axs[3].set_xlabel('Latitude [deg]')
# axs[3].set_ylabel('Pressure [hPa]')

fig.savefig('Figures/comparison_vmr.png', bbox_inches = 'tight')



#-----------------vertical profiles----------------------------------------------------

fig, ax = plt.subplots(1, 1, figsize = [8, 8])
for i in [1, 100, 200]:
    ax.plot(r2vmr[:, i],p_grid,  alpha = 0.8, color = 'c', label = 'r2vmr')
    ax.plot( q2vmr[:, i],p_grid, alpha = 0.8, color = 'r', label = 'q2vmr')
    ax.set_yscale('log')

    ax.set_ylabel('pressure [hPa]')
    ax.set_xlabel('vmr [m3/m3]')   

    ax.set_xscale('log')

    ax.set_ylim(1250, 5)
       
    
    ax2 = ax.twiny()
    color = 'tab:blue'
    ax2.plot(grid_t[:, i], p_grid, 'k', label = 'temperature')
    ax2.set_xlim(180, 260)
ax2.set_xlabel('temp[K]')
ax2.legend( loc = "lower right")
ax.legend(["r2vmr", "q2vmr"])
#    ax2.set_yscale('log')
fig.savefig('Figures/vertical_profile_vmr.png', bbox_inches = 'tight')
