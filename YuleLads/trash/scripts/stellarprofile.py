'''
Plot profiles around stellar particles
Sam Geen, August 2017
'''

from startup import *

import os, glob, time
from adjustText import adjust_text

import Hamu
import pymses
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import linestyles
import sinks, stellars
import skimfunc

import rayprof

import copy

gamma = 1.4 # Ramses default
mH = 1.66e-24
kB = 1.38e-16

def stellarmassinsnap(snap):
    print "Finding stellar mass in snapshot"
    sink = sinks.FindSinks(snap.hamusnap)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    return np.sum(sink.mass)

def maxstellarmassinsnap(snap):
    print "Finding max stellar mass in snapshot"
    stellar = stellars.FindStellar(snap.hamusnap)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    if len(stellar.mass) > 0:
        mmax = np.max(stellar.mass)
    else:
        mmax = 0.0
    return mmax

def runforsim(simname):
    tracks = {}
    sim = Hamu.Simulation(simname)
    boxlen = 0.0
    raylen = 1.0 # pc
    # Make tracks of each star
    snapdict = {}
    print "Setting up tracks"
    for snap in sim.Snapshots():
        snapdict[snap.RawData().iout] = snap
        if boxlen == 0.0:
            boxlen = snap.RawData().info["boxlen"]
        sink = sinks.FindSinks(snap)
        stellar = stellars.FindStellar(snap)
        stellar.FixPosition(sink)
        masses = stellar["mass"]
        for mass in masses:
            slice = stellar.OnlyMass(mass)
            if not mass in tracks:
                tracks[mass] = slice
            else:
                tracks[mass].Append(slice)
    print "Plotting..."
    starid = 0
    for mass, track in tracks.iteritems():
        starid += 1
        for it in range(0,track.Length()):
            centre = [track["x"][it]/boxlen,track["y"][it]/boxlen,track["z"][it]/boxlen]
            length = raylen / boxlen
            snap = snapdict[track["snapid"][it]]
            #print track["snapid"]
            #time.sleep(1)
            rayprof.plot(snap,"rho",centre=centre,length=length,suffix="stellar"+str(starid),
                         simname=simname)
    print "Done!"

def run(simnames):
    for simname in simnames:
        runforsim(simname)

if __name__=="__main__":
    simnames = allsims
    run(simnames)
