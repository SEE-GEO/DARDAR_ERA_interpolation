#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:05:48 2020

@author: inderpreet
"""
import numpy as np
from era2dardar.utils.ellipsoidmodels import ellipsoidmodels
from era2dardar.utils.ellipsoidradii import ellipsoidradii
from typhon import constants
from era2dardar.utils.pos2g import pos2g
from era2dardar.utils.alt2pressure import alt2pres, pres2alt
from scipy.interpolate import interp1d


def alt2pres(altitude):
    '''
    Determine site pressure from altitude.
    Follows z2p_simple from atmlab

    Parameters
    ----------
    altitude : numeric
        Altitude above sea level. [m]

    Returns
    -------
    pressure : numeric
        Atmospheric pressure. [Pa]

    '''

 #   press = 100 * ((44331.514 - altitude) / 11880.516) ** (1 / 0.1902632)
    
    pres = 10. **( 5 - altitude/16e3 )

    return pres

def pres2alt(pressure):
    '''
    Determine altitude from site pressure.
    follows p2z_simple from atmlab

    Parameters
    ----------
    pressure : numeric
        Atmospheric pressure. [Pa]

    Returns
    -------
    altitude : numeric
        Altitude above sea level. [m]

    '''

#    alt = 44331.5 - 4946.62 * pressure ** (0.190263)
    
    alt = 16e3 * ( 5 - np.log10(pressure) )
    

    return alt


def shift2refpoint( p, z, p0, z0 ):
    
    #f = interp1d(np.log(p), z, kind = "linear")
    
    f = interp1d(np.log(p.ravel()), z.ravel(), kind = "linear")    
        
    #f = interpn(np.log(p), z, np.log(p0), method = "linear")
    
    z = z - (f(np.log(p0))  - z0 )
    #z = z - (f  - z0 )

    return z


def z2g(r_geoid, g0, z):

      g = g0 * (r_geoid/(r_geoid+z)) ** 2
      
      return g

def pt2z(p, t, h2o, p0, z0, lat = 45, z_acc = -1):
    """
    ***** Translated from atmlab function pt2z ****
    
    Calculates altitudes fulfilling hydrostatic equilibrium, based on
    vertical profiles of pressure, temperature and water vapour. Pressure
    and altitude of a reference point must be specified.

    Molecular weights and gravitational constants are hard coded and 
    function is only valid for the Earth.

    As the gravitation changes with altitude, an iterative process is
    needed. The accuracy can be controlled by *z_acc*. The calculations
    are repeated until the max change of the altitudes is below *z_acc*. If
    z_acc<0, the calculations are run twice, which should give an accuracy
    better than 1 m.


    Parameters
    ----------
    p : np.array containing pressure values [Pa]
    t : np.array containing temperature [K]
    h2o : np.array or a scalar, Water vapour [VMR]. 
    p0 : Pressure of reference point [Pa]
    z0 : Altitude of reference point [m]
    lat : Latitude [deg]. Default is 45.
    z_acc : Accuracy for z. Default is -1.

    Raises
    ------
    ValueError
        dimensions of p, t should be identical.
        h2o can be scalar but if np.array, it should be of dimension of p

    Returns
    -------
    z, np.array, geometric altitudes fulfilling hydrostatic equilibrium
    """
    
    
    
    
    n = len(p)
    
    if len(t) != n:                                                      
        raise ValueError('The length of *p* and *t* must be identical.')
                       
# Expand *h2o* if necessary                                              
    if  np.isscalar(h2o):
      h2o = np.tile( h2o, n )  
       
    if  len(h2o) != n:                          
            raise ValueError('The length of *h2o* must be 1 or match *p*.')

    #if p0 > p[0]  or  p0 < p[-1]:
    # raise ValueError('Reference point (p0) can not be outside range of *p*.')
      
    if np.any(p0 > p[0])  or  np.any(p0 < p[-1]):
      raise ValueError('Reference point (p0) can not be outside range of *p*.')      
                                                                    

    
    ellipsoid = ellipsoidmodels('wgs84')
    
    
    # make rough estimate of z
    
    z = pres2alt(p)
    z =shift2refpoint(p, z, p0, z0)
    
    # set Earth radius and g at z=0
    
    re = ellipsoidradii(ellipsoid, lat)
    g0 = pos2g(lat, 0)
    
    # gas constants and modelcular weight of dry air and water vapor
    
    r   = constants.R
    md  = 28.966
    mw  = 18.016

    k  = 1-mw/md    # 1 - eps
    rd = 1e3 * r / md  # gas constant for 1 kg of dry air
        
    
    #How to end iterations
    
    if z_acc < 0:
      niter = 2
    else:
      niter = 99
    
    for j in range(niter):
    
          zold = z
        
          g = z2g( re, g0, z )

          for i in range(n-1):

        
                gp  = ( g[i] + g[i+1] ) / 2;
        
                #Calculate average water VMR (= average e/p)
                hm  = (h2o[i]+h2o[i+1]) / 2;
        
                #The virtual temperature (no liquid water)
                tv = (t[i]+t[i+1]) / ( 2 * (1-hm*k) )  # 3.16 in Wallace&Hobbs
        
                #The change in vertical altitude from i to i+1 
                dz = rd * (tv/gp) * np.log( p[i]/p[i+1] )
                z[i+1] = z[i] + dz
            
        
          # Match the altitude of the reference point
          z = shift2refpoint( p, z, p0, z0 );
          if (z_acc >= 0) & (np.max(np.abs(z-zold)) < z_acc):

            break
        
    return z    
