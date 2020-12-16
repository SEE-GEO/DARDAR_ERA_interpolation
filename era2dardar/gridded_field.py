#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 14:19:52 2020


converts gridded fields to ARTS griddedField class
useful to write to xml and use other functions available for
typhon GriddedField class

@author: inderpreet
"""

from typhon.arts.griddedfield import GriddedField1, GriddedField2
from typhon.arts.griddedfield import GriddedField3, GriddedField4


def gridded_field1(field, data, gridnames, grids):
    
    gf             = GriddedField1()
    gf.data        = field
    gf.gridnames   = gridnames
    gf.grids       = grids
    
    return gf

def gridded_field2(field, data, gridnames, grids):
    
    gf             = GriddedField2()
    gf.data        = field
    gf.gridnames   = gridnames
    gf.grids       = grids
    
    return gf

def gridded_field3(field, data, gridnames, grids):
    
    gf             = GriddedField3()
    gf.data        = field
    gf.gridnames   = gridnames
    gf.grids       = grids
    
    return gf

def gridded_field4(field, data, gridnames, grids):
    
    gf             = GriddedField4()
    gf.data        = field
    gf.gridnames   = gridnames
    gf.grids       = grids
    
    return gf