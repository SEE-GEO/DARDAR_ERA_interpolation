#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 15:00:30 2020

interface to generate xml files for multiple DARDARfiles

@author: inderpreet
"""
import os
import glob
import numpy as np
import random

from era2dardar.dardar2atmdata import dardar2atmdata
from era2dardar.DARDAR import DARDARProduct
from era2dardar.atmData import atmdata
from era2dardar.utils.alt2pressure import alt2pres
import typhon.arts.xml as xml
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import shutil
import subprocess
import zipfile
import subprocess
from era2dardar.utils.read_from_zip import read_from_zip 

#from era2dardar.utils.add2zip import add2zip, check_in_zip

def add2zip(zfile, filename):
    
     print (zfile, filename)
    
     
     subprocess.call(["zip", "-j", zfile, filename])


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
    



def run_all_zsurface(p_grid, zipfiles, latlims, zippath):
    
    srtmpath = os.path.expanduser("~/Dendrite/Projects/IWP/GMI/DARDAR_ERA_m65_p65_SRTM/")
    
    for zfile in zipfiles:
        
            dardarfile, N = zip2dardar(zfile)
            
            print (os.path.basename(dardarfile), os.path.basename(zfile))
                
            dardar = DARDARProduct(dardarfile, latlims = latlims, node = N)
            date   = filename2date(dardarfile)
            #dardar.plot_scene()

                
            
            print ("t_0, t_1", dardar.t_0, dardar.t_1)
    
            
            lon          = dardar.longitude 
            lat          = dardar.latitude
            
            lat1         = np.around(lat.min() - 2)
            lat2         = np.around(lat.max() + 2)
            
            lon1         = np.around(lon.min() - 2)
            lon2         = np.around(lon.max() + 2)
    
            if lon1 < -180:
                lon1 = -180.0
            if lon2 > 180:
                lon2 = 180.0    
      
# when encountering prime meridian, download global data
# keeps interpolation simple          
            if lon1 == -180 or lon2 == 180:
                lon1 = -180.
                lon2 =  180.
                
# domain for which ERA5 data is downloaded           
            domain  = [lat1, lat2, lon1, lon2]
            

# instantiate the atmdata class to get only N0star            
            atm        = atmdata(dardar, p_grid, domain = domain)

            z_surface  = atm.z_surface_srtm
            
            parameter  = "z_field"
            zfield    = np.squeeze(read_from_zip(zfile, parameter))
            print(zfield.shape)
            Z          = atm.Z(zfield)
            
            
            print (z_surface.shape)
            
#             year  = date.strftime("%Y")  
#             jday  = date.strftime("%j")
#             hour  = date.strftime("%H")
            
            
#             N0star_path  = os.path.join(zippath, year + jday + hour + N)
            
#             if not os.path.isdir(N0star_path):
#                 os.makedirs(N0star_path)
            zfile_srtm = os.path.join(srtmpath, os.path.basename(zfile))
            
            print(zfile, zfile_srtm)
            
            shutil.copy(zfile, zfile_srtm)

            filename = (os.path.join(srtmpath,  "z_surface" + ".xml"))                
            xml.save(z_surface, filename)
            
            zfieldfile = (os.path.join(srtmpath,  "reflectivities" + ".xml"))                
            xml.save(Z, zfieldfile)
            
            add2zip(zfile_srtm, filename)
            add2zip(zfile_srtm, zfieldfile)
            
            #  remove unzipped folder
            #shutil.rmtree(os.path.dirname(filename))


            
     
            # # remove downloaded ERA files  
#             erafiles = (glob.glob(os.path.join("ERA5/*/", "*" 
#                                                + date.strftime("%4Y") 
#                                                + date.strftime("%2m") 
#                                                + date.strftime("%2d")  
#                                                + date.strftime("%2H") +"*")))
#             for f in erafiles:    
#                 os.remove(f)         
        

if __name__ == "__main__":

    # pressure grid
    p_grid = alt2pres(np.arange(-700, 20000, 250))
    p_grid = np.concatenate([p_grid, np.array([30, 20, 10, 7, 5, 3, 2, 1]) * 100])
    
    # latitudinal extent
    latlims    = [-65, 65]
    
    # year and month of data
#    year = "2009"
#    month = "06"
    
    # dardar inputfile path
#    inpath  = os.path.join(os.path.expanduser("~/Dendrite/SatData/DARDAR"), year, month)
    
    #zipfile path
    zippath = os.path.expanduser("~/Dendrite/Projects/IWP/GMI/DARDAR_ERA_m65_p65_zfield")
    
    zipfiles = glob.glob(os.path.join(zippath, "*.zip"))
    
    #zipfiles = ["/home/inderpreet/Dendrite/Projects/IWP/GMI/DARDAR_ERA_m65_p65_zfield/2010_027_07_A.zip"]
    #zipfiles = ["/home/inderpreet/Dendrite/Projects/IWP/GMI/DARDAR_ERA_m65_p65_zfield/2006_355_07_D_N.zip"]
    
    
    run_all_zsurface(p_grid, zipfiles[700:], latlims,  zippath)
    
    
    

    
        
    #------------------------------------to check data----------------------------    
    # dardar.plot_scene()
    
    # #check IWC values
    # fig, ax = plt.subplots(1,1, figsize = [16, 8])
    # cmap = 'coolwarm'
    # lat_d = dardar.latitude
    # height_d  = dardar.height
    # iwc       = dardar.iwc
    # p         = alt2pres(height_d) * 0.01
    # im = ax.pcolormesh(lat_d, p, iwc.T, cmap=cmap, norm=colors.LogNorm(),
    #                                       vmin = 0.000001, vmax = 0.01)
    # ax.set_ylim(ax.get_ylim()[::-1])
    # # #ax.set_yscale('log')
    # fig.colorbar(im, ax=ax, label = 't [K]' , extend = 'both')
    # #
    # # get atmdata on DARDAR grid
    
    # atm = atmdata(dardar, p_grid)
    
    # #lsm = atm.sea_ice_cover
    # #plt.plot(lat_d, lsm)
    # z = atm_fields["z_field"]
    # fig, ax = plt.subplots(1,1, figsize = [8, 8])
    # ax.plot(np.log(p_grid * 0.01), z[:, 1000, 0]/1000)
    # ax.plot(np.log(p_grid * 0.01), z[:, 2000, 0]/1000)
    
    # ax.plot(np.log(p_grid * 0.01), z[:, 3, 0]/1000)
    
    # fig, ax = plt.subplots(1, 1, figsize = [8,8])

    # ax.plot(lat, lsm)
    # ax.plot(lat[mask], lsm[mask], )
    # ax.plot(lat[~mask], lsm[~mask], )
