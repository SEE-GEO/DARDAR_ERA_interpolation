#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 08:42:59 2020

@author: inderpreet
"""


def scale_vmr(vmr_n2, vmr_h20):
    """
    

    Parameters
    ----------
    vmr_n2  (ndarray) : array containing standard vmr values for N2/O2
    vmr_h20 (ndarray) : array containing vmr for H2O, for the same 
    grid as vmr_n2

    Returns
    -------
    vmr_n2 (ndarray)  : scaled vmr for N2/O2

    """
    vmr_n2  = vmr_n2 * (1 - vmr_h20)
    
    
    return vmr_n2