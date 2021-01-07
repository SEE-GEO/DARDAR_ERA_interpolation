#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 10:44:42 2020

convert pressure to altitude and vice-versa

"""
import numpy as np

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
