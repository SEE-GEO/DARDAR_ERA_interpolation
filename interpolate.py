#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

A general purpose 3d interpolator function for 
data on a regular grid

Created on Wed Dec  2 11:09:54 2020

@author: inderpreet
"""

from scipy.interpolate import RegularGridInterpolator

def interpolate(level, lat, lon, A):
        """
        interface to scipy 3D linear interpolator
        All dimensions should be in ascending order

        Parameters
        ----------
        level : np.array of dimension P, gives the z coordinate
        lat : np.array of dimension M
        lon : np.array of dimension N
        A : field values to be interpolated from, dimension PxNxM

        Returns
        -------
        interpolator function

        """

        return (RegularGridInterpolator((level, lat, lon), A, 
                                        bounds_error = True))