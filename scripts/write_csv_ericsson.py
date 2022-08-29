#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 11:04:58 2021

@author: inderpreet
"""
import os
import numpy as np
import pandas as pd
from typhon.physics.thermodynamics import e_eq_water_mk
from era2dardar.utils.refractive_index import refractive_index
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from era2dardar.utils.alt2pressure import pres2alt
import pickle
plt.rcParams.update({'font.size': 16})

#%%

#----------------------------------------------------------------
# reading Arlanda data and computing refractive index
df = pd.read_csv('arlanda_data.csv', skipinitialspace = True)


p  = df["PRES"]
T  = df["TMP"]
Td = df["DPT"]
H  = df["HGT"] 

e_w = e_eq_water_mk(Td)

n = refractive_index(p, e_w, T)
n_ppm = (n -1) * 1000000

#------------------------------------------------------------------


#t0 = datetime(2021, 7, 23, 23, 00, 00 )
#t1 = datetime(2021, 7, 29, 23, 00, 00)

t0 = datetime(2020, 9, 10, 23, 00, 00 )
t1 = datetime(2020, 9, 12, 00, 00, 00)
    
# location A
#lat = 59.86
#lon = 17.76

# location B
#lat = 59.50
#lon = 18.06    
key = "A"
t2 = t1 
while t0 < t1:
        
    t0 = t0 + timedelta(minutes = 60)   
    with open("arlanda_era5_" + str(t0.day) + str(t0.hour) + "_" + key + ".pickle", 'rb') as f:
        grids = pickle.load(f)
        p_era = pickle.load(f)
        f.close()
        
    print ("arlanda_era5_" + str(t0.day) + str(t0.hour) + "_" + key + ".pickle") 
    T_era = grids["t"]
    r_era = grids["r"]
    
    e_w_era = r_era * e_eq_water_mk(T_era)/100
    
    n_era = refractive_index(p_era * 100, e_w_era.ravel(), T_era.ravel())
    
    n_era_ppm = (n_era -1) * 1000000
    
    era_h = pres2alt(p_era * 100)
    
    loc_data = {}
    
    loc_data["height"] = era_h[:-1]
    loc_data["pressure"] = p_era[:-1]
    loc_data["relative_humidity"] = r_era[:-1, 0]
    loc_data["temperature"] = T_era[:-1, 0]
    
    loc_data["refractive_index"]  = n_era_ppm[:-1]

    df = pd.DataFrame.from_dict(loc_data)
    
    local_time = t0 - timedelta(minutes = 60)
    outfile = "loc_" + key + "_" + local_time.strftime("%Y%m%d:%H")
    print (outfile)
    df.to_csv("refractive_index_data_ERA5/" + outfile + ".csv")
    
    




 
