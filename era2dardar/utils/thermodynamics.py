#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 12:47:25 2020

@author: inderpreet
"""
import typhon.physics.thermodynamics as thermodynamics 
import typhon.physics.atmosphere as atmosphere
import numpy as np
from era2dardar.utils.alt2pressure import pres2alt 
import scipy
from era2dardar.utils.integrate import integrate


def mixr2massconc( mixr, pres, temp):
    """
    convert mass mixing ratio to mass concentratio 

    Parameters
    ----------
    mixr : mass mixing ratio [kg/kg]
    pres : corresponding pressure value [Pa]
    temp : corresponding temperature value [K]

    Returns
    -------
    massconc : mass concentration  [kg/m3]

    """
    rho =  airdensity(pres, temp ) # kg / m3
    
    #= mass concentration
    massconc = rho * mixr  # kg / m3
    
    return massconc


def airdensity(p,t):
    """
    calculate density of air at given pressure and temperature

    Parameters
    ----------
    p : pressure [Pa]
    t : temperature [K]

    Returns
    -------
    rho : air density [kg/m3]

    """
    Rsp = 287.058;  # Jkg^{-1}K^{-1}
    
    rho = p / ( Rsp * t ) 
    
    return rho

def vmr2mixr(vmr, M_w):
    """
    convert volume mixing ratio to mass mixing ratio
    
    Parameters
    ----------
    vmr : np.array, volume mixing ratio [m3/m3]
    M_w : modelcular mass of species [kg/mol]
    M_d : modelcular mass of dry air [kg/mol] 

    Returns
    -------
    mass mixing ratio [Kg/Kg]

    """
    
    M_d = 28.9661e-3
    
    mixr = vmr * M_w/M_d
    
    return mixr


def mixr2vmr( mixr, M_w ):
    """
    convert mixing ratio to volume mixing ratio

    Parameters
    ----------
    mixr : np.array, mass mixing ratio [Kg/Kg]
    M_w : modelcular mass of species [kg/mol]
    M_d : modelcular mass of dry air [kg/mol] 

    Returns
    -------
    vmr : volume mixing ratio [m3/m3]

    """
    

    # molecular mass of dry air
    M_d = 28.9661e-3;
    
    vmr = mixr * ( M_d / M_w )
    return vmr


def rh2vmr(grid_p, grid_t, grid_r):
    """
    conversion of relative humidty values to volume mixing ratio (vmr)

    Parameters
    ----------
    p_grid : [n x m] np.array, pressure grid defining n pressure levels at m locations [Pa] 
    grid_t : [n x m] np.array, temp gridded to DARDAR lat locations and levels defined 
    by p_grid [K]
    grid_r : [n x m] np.array, relative humidity gridded to DARDAR lat locations and levels defined 
    by p_grid [%]

    Returns
    
    np.array, volume mixing ratio [m3/m3]

    """

    
    t_liq       = grid_t >= 273.0
    t_ice       = grid_t <= 240.0
    t_mix       = ( grid_t > 240.0) & (grid_t < 273.0 )
    
    grid_r      = grid_r/ 100
    
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
    
    return r2vmr

def massconc2tcwv(mc, p_grid, grid_sp):
    """
    integrate the mass concentration to get the TCWV over each grid point

    Parameters
    ----------
    mc : np.array [mxn] mass concentration of species
    p_grid : np.array, m pressure leveles [Pa]
    grid_sp : np.array [n], surface pressure [Pa]

    Returns
    -------
    tcwv : TCWV  [kg/m2]

    """

    N        = grid_sp.shape[0]    
    tcwv     = np.zeros(N)
    for i in range(N):
        inds = np.where(p_grid < grid_sp[i])[0]
        pres = np.append(p_grid[inds], grid_sp[i])
        
        A = mc[inds, i]
        A = np.append(A, mc[-1, i])
        
        tcwv[i] = integrate(A, pres2alt(pres))
        
    return tcwv 

