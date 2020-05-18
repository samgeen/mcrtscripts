'''
Find stellar tempertaure from mass
Sam Geen, February 2019
'''

import os, sys
sys.path.append("../StellarSources")
sys.path.append("../StellarSources/Fortran/f2py")

import numpy as np
import scipy.interpolate

import cPickle as pik

from Input import readgeneva, derivedgeneva, readstarburst, readsbwinds

import singlestar, makeimf
from IMFUtils import sampleimf

from consts import *

# Verbose mode?
verbose = True

# Stellar track stuff
windtracks = None
trackloc = "../StellarSources/"
starmasses = np.arange(5,121,5)
starmetals = [0.002,0.014]

# Set up single star module
# startableloc = "/home/samgeen/Programming/Astro/StellarSources/Compressed/singlestar_z0.014"
# print "Reading star tables in", startableloc
# singlestar.star_setup(startableloc)  

time = 1e6 # yr
f = open("MtoTeff.txt","w")
for metal in starmetals:
    print "Z =", metal
    f.write("---\n")
    f.write("Z ="+str(metal)+"\n")
    f.write("M / Msun | Teff / K\n")
    for mass in starmasses:
        print mass,"..."
        tracks = readgeneva.Tracks()
        track = tracks.FindTrack(mass,metal,rotating=True)
        spectrack = readstarburst.Track(mass,metal,trackloc)
        spectimes = spectrack.Times() # in yr
        QHs = np.array(spectrack.NPhotonsVsTime(13.6, 1000.0))
        Tes = 10.0**track["logTe"]
        ts = np.abs(track["Time"] - time) # in yr
        sts = np.abs(spectimes - time) # in yr
        Te = Tes[ts == ts.min()][0]
        QH = QHs[sts == sts.min()][0]
        f.write(str(mass)+" "+str(Te)+" "+str(QH)+"\n")
    f.write("---\n")
f.close()