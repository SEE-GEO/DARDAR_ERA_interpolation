#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 12:50:51 2020

@author: inderpreet
"""
import numpy as np

def lon_180_to_360(lon):
    return lon % 360
        

def lon_360_to_180(lon):
    return(lon + 180) % 360 - 180