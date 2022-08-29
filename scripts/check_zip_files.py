#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 13:23:56 2021

@author: inderpreet
"""


import os
import glob
from era2dardar.utils.add2zip import  check_in_zip


def check_zip_file(zfile, infiles):
    
    for file in infiles:
        
        flag = check_in_zip(zfile, file)
        if flag == False:
            print (flag)
            print (zfile)
            return
        
        
        
        
infiles =               ["time.xml",              
                        "lon_grid.xml",
                        "lat_grid.xml",
                        "p_grid.xml" ,
                        "t_field.xml" ,
                        "z_field.xml",
                        "lwc.xml" ,
                        "skt.xml" ,
                        "t2m.xml",
                        "wind_speed.xml" ,
                        "wind_direction.xml" ,
                        "iwc.xml" ,
                        "z_surface.xml" ,
                        "abs_species.xml" ,                   
                        "vmr_field.xml", 
                        "lsm.xml"  ,
                        "sea_ice_cover.xml" ,
                        "snow_depth.xml" ,
#                        "reflectivities.xml" ]
                        "N0star.xml"]      






zfiles = glob.glob("/home/inderpreet/Dendrite/Projects/IWP/IceCube/DARDAR_ERA_m30_p30_N0star/*.zip")
#zfiles = glob.glob("/home/inderpreet/Dendrite/Projects/IWP/GMI/DARDAR_ERA_m65_p65/*.zip")

for zfile in zfiles:
 #   print (zfile)
    
    check_zip_file(zfile, infiles)