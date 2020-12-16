#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 21:16:24 2020

@author: inderpreet
"""


import datetime as dt

def seconds2datetime(n): 
    """
    We use datetime.timedelta 
    to represent seconds into 
    hours, minutes and seconds format

    Parameters
    ----------
    n :  scalar, seconds

    Returns
    -------
    datetime object

    """
    time = str(dt.timedelta(seconds = n)) 
    datetime_obj = dt.datetime.strptime(time, '%H:%M:%S.%f')
    
    return datetime_obj
    