'''
Code to run on script startup for this project
Includes imports that can be used by other scripts
Sam Geen, June 2016
'''

import sys, os, glob
import numpy as np
# Import from the main scripts folder
sys.path.append("../../scripts")

import customplot
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import hydrofuncs
import snaptime
import shade4Myr

import linestyles

import Hamu
import pymses

Hamu.Workspace("HIISFE")

# Physical conversions
X = 0.76
kB = 1.38062e-16
pcincm = 3.086e18
Msuning = 1.9891e33
mHing = 1.66e-24
Myrins = 3.1556926e13

# Some global definitions
Msolar = "M$_{\odot}$"

# Simulation names
allsims = ["L-RT","L-NRT",
           "M-RT","M-NRT",
           "S-RT","S-NRT",
           "XS-RT","XS-NRT"]

# Useful functions
# N_H to/from A_k (from Lombardi+ 2010)
akconst = 1.75e22
def AktoNH(Ak):
    return akconst*Ak

def NHtoAk(NH):
    return NH/akconst

def ColDenstoNH(dens):
    return dens*X/mHing

def NHtoColDens(NH):
    return NH*mHing/X

def ColDenstoAk(dens):
    return NHtoAk(ColDenstoNH(dens))

def AktoColDens(Ak):
    return NHtoColDens(AktoNH(Ak))

# Conversion constant from Kirk et al 2013
muH2 = 2.83
def NH2toColDens(NH2):
    return NH2*mHing*muH2

def ColDenstoNH2(dens):
    return dens/(mHing*muH2)

def NH2toNH(NH2):
    return NH2*X*muH2

def NHtoNH2(NH):
    return NH/(X*muH2)

def MakeDirs(folder):
    try:
        print "Making directory", folder
        os.makedirs(folder)
    except:
        print "Not making directory", folder, "- it already exists!"
        pass

print "Imported various project-global modules"
