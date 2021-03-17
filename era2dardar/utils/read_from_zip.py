#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 21:14:35 2021

@author: inderpreet
"""

from typhon.arts import xml
import zipfile
import numpy as np
def read_from_zip(zfile, parameter):

    with zipfile.ZipFile(zfile, 'r') as zip_ref:
        zip_ref.extractall("temp/")
    
    var = np.squeeze(xml.load("temp/" + parameter + ".xml"))  
    return var