#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:28:45 2021

@author: inderpreet
"""


def refractive_index(p, e, T):
    """
    The (non-dispersive) refractive index n is calculated using Eq 8 of 
    Buehler et al 2005, JQRST


    Parameters
    ----------
    p : total air pressure [Pa]
    e : water vapour pressure [Pa]
    T : temperature [K]

    Returns
    -------
    n : (non-dispersive) refractive index

    """
    n = 1 + 77.593 * 1e-8 * (p - e)/T + e * (72 * 1e-8/T  +  3.754 * 1e-3/T**2)
    
    return n
 