#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 12:43:37 2021

@author: inderpreet
"""
import zipfile
import os

def add2zip(zfile, filename):
    """
    Adds a new file to existing zipfolder.

    Parameters
    ----------
    zfile : path to existing zip folder
    filename : path to the file to be added

    Raises
    ------
    Exception
        if file already exists in the zip file

    Returns
    -------
    None.

    """


    if check_in_zip(zfile, os.path.basename(filename)):
        raise Exception("File already exists, cannot overwrite")
        
    with zipfile.ZipFile(zfile, 'a') as zf:

        source_path = filename
        destination = os.path.basename(filename)
        zf.write(source_path, destination)
        zf.close()
        
    
def check_in_zip(zfile, filename):
    """
    checks if a file exists in the zipfile

    Parameters
    ----------
    zipfile : string containing name of zipfile with path
    filename : string, file to be checked.

    Returns
    -------
    bool
        True: if file exists
        False : otherwise

    """
    zf = zipfile.ZipFile( zfile, "r", zipfile.ZIP_DEFLATED)

    if filename in zf.namelist():
        zf.close()  
        return True
    else:
        zf.close()  
        return False