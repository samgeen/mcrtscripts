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

import linestyles

import HamuLite as Hamu
import pymses

import sinks
import stellars
import timefuncs

#Hamu.Workspace("HIISFE")

# Physical conversions
X = 0.76
kB = 1.38062e-16
pcincm = 3.086e18
Msuning = 1.9891e33
mHing = 1.66e-24
Myrins = 3.1556926e13

# Simulation units (to prevent need to dig into code)
unit_l      =  0.308000000000000E+19
unit_d      =  0.232474340000000E-23 
unit_t      =  0.253950790000000E+16

# Simulation global properties
tffcloud_Myr = 4.21875 # tff of cloud as a whole in Myr
tffcloud_code = tffcloud_Myr / (unit_t/Myrins)

# Some misc definitions
Msolar = "M$_{\odot}$"

# Plot folder
plotfolder = "../plots/"

# Simulation names
imfsims = ["IMF"+str(i).zfill(2) for i in range(1,14)]
icsims  = ["IC"+str(i).zfill(2) for i in range(1,14)]
allsims = imfsims + icsims

# Simulation locations
#mainsimfolder = "/home/hd/hd_hd/hd_mp149/MCRT/runs/"
mainsimfolder = "/home/sgeen/MC_RT/runs_anais/"
imfsimfolder = mainsimfolder+"/55_Lcloud_yulelads/"
icsimfolder = mainsimfolder+"/57_Lcloud_ictest/"

# Set up simulation arrays
def _MakeSim(name):
    simexists = True
    try:
        sim = Hamu.Simulation(name,silent=True)
    except KeyError:
        simexists = False
    if not simexists:
        # Identify simulation folder
        if "IMF" in name:
            folder = imfsimfolder
        elif "IC" in name:
            folder = icsimfolder
        else:
            print "Not IMF or IC in sim name", name
            raise ValueError
        folder += name[-2:]+"_justuv/"
        if not os.path.exists(folder):
            print "No folder", folder
            raise OSError
        # Return label
        # TODO! FOR NOW, JUST RETURN NAME IMF/IC+NUM
        label = name
        sim = Hamu.MakeSimulation(name,folder,label)
    return sim

hamusims = {}
for simname in allsims:
    hamusims[simname] = _MakeSim(simname)


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