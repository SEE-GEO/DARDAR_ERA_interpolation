#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 20:54:46 2021

@author: inderpreet
"""
import os
import glob
from datetime import datetime

def filename2date(filename):
        filename = os.path.basename(filename)
        filename = filename.split("_")[2]
        pattern = "%Y%j%H%M%S"
        return datetime.strptime(filename, pattern)
    
def find_dardar(date):
    """
    get dardar filename from the zipfile name

    Parameters
    ----------
    date : datetime object containing timestamp of DARDAR file
    Returns
    -------
    dardarfile : string containing DARDAR filename

    """
    
    inpath = os.path.expanduser("~/Dendrite/SatData/DARDAR")
    
    year  = date.strftime("%Y")
    month = date.strftime("%m")    
    jday  = date.strftime("%j")
    day   = date.strftime("%d")
    hour  = date.strftime("%H")
    
    dardarfile = (glob.glob(os.path.join(inpath, year, month, day , 
                                         "*" + year + jday + hour + "*"  )))
    
    return dardarfile[0]
    
def zip2dardar(zfile):
    
    basefile   = os.path.basename(zfile)
    filename   = basefile[:11]
    pattern    = "%Y_%j_%H"
    
    date       = datetime.strptime(filename, pattern)

    if "A" in basefile:
        N = "A"
    if "D_S" in basefile:
        N = "D_S"
    if "D_N" in basefile:
        N = "D_N"        
    
    dardarfile = find_dardar(date)
    
    return dardarfile,  N