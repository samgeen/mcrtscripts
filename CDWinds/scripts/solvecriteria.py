'''
Calculate point at which stars should overrun their photoionisation bubble with winds/radiation pressure
Sam Geen, December 2017 (lol, try August 2018)
'''

import sys

import customplot

import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate
import ionisedtemperatures

from consts import *

#import stars
#from stars import _Star # needed to make caching work

from collections import OrderedDict

# Stellar track stuff
windtracks = None
trackloc = "../StellarSources/"

Msolar = "M$_{\odot}$"

clabelfmt = '%1.1f'

# Colour maps for the various variables
cmaps = {}
cmaps["Rstall"] = "BrBG"
cmaps["Cb"] = "PiYG_r"
cmaps["Cw"] = "RdBu_r"
cmaps["Crp"] = "PuOr_r"

def SoundSpeed(cloud,star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),cloud.metal)
    # Find ci
    # 2.0 factor is because of free electrons
    ci = np.sqrt(2.0 * gamma * kB * Ti * X / mH)
    return ci

def SoundSpeedFromT(temp):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = temp
    # Find ci
    # 2.0 factor is because of free electrons
    ci = np.sqrt(2.0 * gamma * kB * Ti * X / mH)
    return ci

def FindalphaB(cloud,star,t):
    #HACK - ignore change in Ti due to ionised gas density
    Ti = ionisedtemperatures.FindTemperature(star.Teff(t),cloud.metal)
    return alpha_B_HII(Ti)

