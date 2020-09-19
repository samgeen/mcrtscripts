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
import singlestar

from collections import OrderedDict

verbose = True

if verbose:
    print("Imports done...")
#Hamu.Workspace("HIISFE")

# Which machine are we on?
myMachine = os.environ["MYMACHINE"]

# Physical conversions
X = 0.76
mH = 1.6735326381e-24
kB = 1.38062e-16
G = 6.674e-8
gamma = 1.4 # RAMSES hard-coded
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

# Set up simulations
if verbose:
    print("Making sims...")

# Simulation locations
#mainsimfolder = "/home/hd/hd_hd/hd_mp149/MCRT/runs/"
#mainsimfolder = "/home/sgeen/MC_RT/runs_anais/"
mainsimfolder = "/greenwhale/LEGO/"
Hamu.CACHEPATH = "/greenwhale/samgeen/cache/"
allsims = {"128_LEGO":"128_LEGO_HM6Z002SNLT",}

if myMachine == "CARTESIUS":
    print("Setting up for Cartesius")
    Hamu.CACHEPATH = "/home/samgeen/amun/cache/"
    mainsimfolder="/home/samgeen/amun/runs/"
    allsims = {"UVWIND_30":"69_AMUN_onestar/13_uv+winds_30",
               "UVWIND_60":"69_AMUN_onestar/11_uv+winds_60",
               "UVWIND_120":"69_AMUN_onestar/03_uv+winds_120"}

hamusims = OrderedDict()
for simname, folder in allsims.items():
    label = simname
    sim = Hamu.MakeSimulation(simname,mainsimfolder+folder,label)
    hamusims[simname] = sim

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
            print("Not IMF or IC in sim name", name)
            raise ValueError
        folder += name[-2:]+"_justuv/"
        if not os.path.exists(folder):
            print("No folder", folder)
            raise OSError
        # Return label
        # TODO! FOR NOW, JUST RETURN NAME IMF/IC+NUM
        label = name
        sim = Hamu.MakeSimulation(name,folder,label)
    return sim


# Set up single star module
if verbose:
    print("Singlestar Setup...")
startableloc = "/home/samgeen/StellarSources/Compressed/singlestar_z0.002"
singlestar.star_setup(startableloc)

# Useful functions
MW = False
SMC = True
if MW:
    # N_H to/from A_k (from Lombardi+ 2010)
    akconst = 1.75e22
    # N_H to/from A_v (from Savage & Mathis 1979)
    avconst = 1.87e21
elif SMC:
    # Dobashi et al. 2009 (last two paragraphs of section 5.2)
    avconst = 1.0/1.5e-22
    # Hack - just use 10x the avconst for Ak
    akconst = 10.0*avconst

def AktoNH(Ak):
    return akconst*Ak

def NHtoAk(NH):
    return NH/akconst

def AvtoNH(Av):
    return avconst*Av

def NHtoAv(NH):
    return NH/avconst

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
        print("Making directory", folder)
        os.makedirs(folder)
    except:
        print("Not making directory", folder, "- it already exists!")
        pass

def alpha_B_HII(T):
    # HII recombination rate
    # input  : T in K 
    # output : HII recombination rate (in cm3 / s)
    l = 315614./T
    a = 2.753e-14 * l**1.5 / (1. + (l/2.74)**0.407)**2.242
    return a   

print("Imported various project-global modules")
