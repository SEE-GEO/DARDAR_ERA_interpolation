#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 12:33:03 2020

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
from era2dardar.utils.mix2massconc  import mixr2massconc
from integrate_vmr import integrate_vmr
plt.rcParams.update({'font.size': 18})




#p_grid = None
#print(p_grid)


# DARDAR file
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/11/06/DARDAR-CLOUD_v2.1.1_2015310114257_50670.hdf")
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071190249_47194.hdf")
file = os.path.expanduser("~/Dendrite/SatData/DARDAR/2015/03/12/DARDAR-CLOUD_v2.1.1_2015071104824_47189.hdf")

filename = file

# create DARDAR instance
dardar = DARDARProduct(filename, latlims = [-5 , 5])
  
# time stamp of DARDAR data 
t_0 = dardar.filename2date()
t_1 = t_0 + timedelta(minutes = 30)


var         = "temperature"
ERA_t       = ERA5p(t_0, t_1, var)


# pressure grid
z_grid      = np.arange(-100, 25000, 250) 
p_grid      = ERA_t.era["level"].data[5:]
#p_grid      = alt2pres(z_grid) * 0.01

# temperature
var         = "temperature"
ERA_t       = ERA5p(t_0, t_1, var)
grid_t      = ERA_t.interpolate(dardar, p_grid)


z_grid      = pres2alt(p_grid * 100)
# pressure grid
grid_p      = np.tile(p_grid, (grid_t.shape[1], 1))
grid_p      = grid_p.T * 100


var         = "specific_humidity"
ERA_q       = ERA5p(t_0, t_1, var)
grid_q      = ERA_q.interpolate(dardar, p_grid)
q2vmr       = thermodynamics.specific_humidity2vmr(grid_q)


var         = "relative_humidity"
ERA_r       = ERA5p(t_0, t_1, var)
grid_r      = ERA_r.interpolate(dardar, p_grid) /100

t_liq       = grid_t >= 273.0
t_ice       = grid_t <= 240.0
t_mix       = ( grid_t > 240.0) & (grid_t < 273.0 )


r2vmr         = np.zeros(grid_r.shape)
r2vmr[t_liq]  = (atmosphere.relative_humidity2vmr(grid_r[t_liq],
                                                        grid_p[t_liq], 
                                                        grid_t[t_liq], 
                                                        thermodynamics.e_eq_water_mk))
r2vmr[t_ice]  = (atmosphere.relative_humidity2vmr(grid_r[t_ice],
                                                        grid_p[t_ice], 
                                                        grid_t[t_ice], 
                                                        thermodynamics.e_eq_ice_mk))
r2vmr[t_mix]  = (atmosphere.relative_humidity2vmr(grid_r[t_mix],
                                                        grid_p[t_mix], 
                                                        grid_t[t_mix], 
                                                        thermodynamics.e_eq_mixed_mk))

#grid_lwc    = grid_lwc * rho  # convert lwc units to mass concentration                                                        thermodynamics.e_eq_mixed_mk))
var             = "total_column_water_vapour"
ERA_tcwv        = ERA5s(t_0, t_1, var)
grid_tcwv       = ERA_tcwv.interpolate(dardar)


#--------------------convert vmr to tcwv------------------------------------------
N_r = mixr2massconc( r2vmr, grid_p, grid_t)
N_q = mixr2massconc( q2vmr, grid_p, grid_t)

tcwv_q = integrate_vmr(N_r[::-1, :], z_grid[::-1])
tcwv_r = integrate_vmr(N_q[::-1, :], z_grid[::-1])

d_r = (grid_tcwv - tcwv_r)/grid_tcwv
d_q = (grid_tcwv - tcwv_q)/grid_tcwv

fig, ax = plt.subplots(1, 1)
ax.plot(d_r)
ax.plot(d_q)
#----------------------plotting data below ------------------------------------------

lat_d    = dardar.get_data('latitude')
fig, axs = plt.subplots(4, 1, figsize = [20, 25])

fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2, hspace = 0.5)

#for var, ax in zip([r2vmr, q2vmr], axs.ravel()):
for var, ax in zip([r2vmr, q2vmr], axs.ravel()):
    im = (ax.pcolormesh(lat_d, p_grid, var, shading = 'auto', cmap = 'tab20c',
                        vmin = 0.000001, vmax = 0.025, 
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
#d[t_mix] = np.nan
#d[t_liq] = np.nan
im = (axs[2].pcolormesh(lat_d, p_grid, d, shading = 'auto', cmap = 'coolwarm',
                        vmin = -20, vmax = 20))

fig.colorbar(im, ax=axs[2], extend = 'both')
axs[2].set_xlabel('Latitude [deg]')
axs[2].set_ylabel('Pressure [hPa]')
axs[2].set_yscale('log')
axs[2].set_ylim(axs[2].get_ylim()[::-1])

# t = grid_t.copy()
# t[t_mix] = np.nan

im = (axs[3].pcolormesh(lat_d, p_grid, grid_t, shading = 'auto', cmap = 'coolwarm', ))
axs[3].set_ylim(axs[3].get_ylim()[::-1])
fig.colorbar(im, ax=axs[3], extend = 'both')
axs[3].set_yscale('log')
axs[3].set_title('temperature')
axs[3].set_xlabel('Latitude [deg]')
axs[3].set_ylabel('Pressure [hPa]')

# im = (axs[3].pcolormesh(lat_d, alt2pres(dardar.get_data('height'))* 0.01, 
#                   dardar.get_data('iwc').T, shading = 'auto', cmap = 'viridis', 
#                   norm=colors.LogNorm()))
# axs[3].set_ylim(axs[3].get_ylim()[::-1])
# fig.colorbar(im, ax=axs[3], extend = 'both')
# axs[3].set_yscale('log')
# axs[3].set_title('IWC')
fig.savefig('Figures/comparison_vmr.png', bbox_inches = 'tight')


fig, ax = plt.subplots(1, 1, figsize = [8, 8])
for i in [500]:
    ax.plot(r2vmr[:, i],p_grid,  alpha = 0.8, color = 'c', label = 'r2vmr')
    ax.plot( q2vmr[:, i],p_grid, alpha = 0.8, color = 'r', label = 'q2vmr')
    ax.set_yscale('log')

    ax.set_ylabel('pressure [hPa]')
    ax.set_xlabel('vmr')   

    ax.set_xscale('log')
    ax.legend(['r2vmr', 'q2vmr'])
    ax.set_ylim(1250, 0.7)
fig.savefig('Figures/vertical_profile_vmr.png', bbox_inches = 'tight')
