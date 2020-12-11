#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 21:22:37 2020

@author: inderpreet
"""

import scipy

def integrate(values, points):    
    
    if list(points) == sorted(points,reverse=True):
        values = values[::-1]
        points = points[::-1]
        
    tcwv = scipy.integrate.trapz(values, points)
    
    return tcwv