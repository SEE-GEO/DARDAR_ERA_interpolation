#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

it is simply an interface to scipy.RegularGridInterpolator
for Interpolation on a regular grid in arbitrary dimensions

Created on Wed Dec  2 11:09:54 2020

@author: inderpreet
"""

from scipy.interpolate import RegularGridInterpolator, interpn

def interpolator(points, A, method):
        """
        interface to scipy nD linear interpolator
        All dimensions should be in ascending order

        Parameters
        ----------
        points: tuple of ndarray of float, with shapes (m1, ), …, (mn, )
        A : field values to be interpolated from, dimension PxNxM
        
        method : "linear", "nearest"

        Returns
        -------
        interpolator function

        """
        return (RegularGridInterpolator(points, A, 
                                        bounds_error = False, method = method, 
                                        fill_value = None))
