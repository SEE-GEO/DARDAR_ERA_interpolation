#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 12:47:25 2020

@author: inderpreet
"""



def mixr2massconc( mixr, pres, temp):

#= density of dry air
    rho =  airdensity(pres, temp ) # kg / m3
    
    #= mass concentration
    massconc = rho * mixr  # kg / m3
    
    return massconc


def airdensity(p,t):

    Rsp = 287.058;  # Jkg^{-1}K^{-1}
    
    rho = p / ( Rsp * t ) 
    
    return rho

