#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 12:34:54 2020

@author: inderpreet
"""
import scipy
import numpy as np

def integrate_vmr(vmr, p_grid):
    
    tcwv = scipy.integrate.trapz(vmr, p_grid)
    
    return tcwv
    