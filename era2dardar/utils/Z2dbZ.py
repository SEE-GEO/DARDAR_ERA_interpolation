#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 11:46:34 2020

@author: inderpreet
"""
import numpy as np

def Z2dbZ(Z):
    """
    convert Z to dBZ

    Parameters
    ----------
    Z : np.array, reflectivities [mm6 m-3]

    Returns
    -------
    dbZ : np.array, decibels of Z
    """
    
    return 10 * np.log10(Z)


def dbZ2Z(dbZ):
    """
    convert dBZ to Z

    Parameters
    ----------
    dbZ : np.array, reflectivities [dbZ]

    Returns
    -------
    Z : np.array, reflectivities [mm6 m-3]
    """
    
    
    return 10 ** (dbZ/10)