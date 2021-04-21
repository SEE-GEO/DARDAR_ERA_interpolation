#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:54:20 2021

@author: inderpreet
"""
import os
import glob
from datetime import datetime

def match_dardar_cloudsat(dardarfiles, cpath):
    """
    matches dardar and cloudsat 2bgeoprof files

    Parameters
    ----------
    dardarfiles : list of dardarfiles
    cpath : path of cloudsat files
    Returns
    -------
    matched_dfiles :list of dardarfiles matching cloudsat files
    cfiles : list of matching cloudsat files

    """
    
    matched_dfiles = []
    cfiles = []
    for dfile in dardarfiles:
        bname    = os.path.basename(dfile)
        granule  = bname.split("_")[3][:5]
        datestr  = bname.split("_")[2]
        pattern  = "%Y%j%H%M%S"
        date     = datetime.strptime(datestr, pattern)
        yday     = date.strftime("%j")
        year     = date.year
        cfile    = glob.glob(os.path.join(cpath, str(year), yday, '*' + granule + '*'))
         
        if len(cfile) != 0:
            cfiles.append(cfile[0])
            matched_dfiles.append(dfile)  
            
    return matched_dfiles, cfiles     