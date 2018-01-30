'''
Generic plotting of runs as a time series
Sam Geen, February 2014
'''
                                                                                                                                                  
import matplotlib as mpl
mpl.use('Agg')

import os
import Hamu
import numpy as np
import matplotlib.pyplot as plt

import matzner2002

import pymses
from pymses.filters import CellsToPoints
from pymses.utils import constants as C
colconv = mpl.colors.ColorConverter()

from collections import OrderedDict

sims = OrderedDict()

import rayprof

def radialstats(snap):
    '''
    Stats on radius of shockwave with angle
    '''
    nprofs, nray = rayprof.nprofs, rayprof.nray
    r, p = rayprof.makeprofsHamu(snap.hamusnap,"xHII")
    radii = r[0:nray]
    profs = np.reshape(p,(nprofs, nray)).T # No idea
    ilim = 0.9
    if profs.max() < ilim:
        return np.zeros((5))
    rmaxes = np.zeros((nprofs))
    # Too tired to think of a way to do this without a for loop
    for iprof in range(0,nprofs):
        lim = np.max(np.where(profs[0:nray,iprof] > ilim))
        rmaxes[iprof] = radii[lim]
        print radii.shape, lim
    percs = np.percentile(rmaxes,[0,25,50,75,100])
    return percs

radialstatsHamu = Hamu.Algorithm(radialstats)

def funcovertime(sim, func):
    times = sim.Times()
    utime = sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)
    times = times * utime - 1.25
    results = np.zeros((len(times),5))
    funcHamu = Hamu.Algorithm(func)
    for itime, snap in zip(range(0,len(times)),sim.Snapshots()):
        results[itime,0:5] = funcHamu(snap)
    results = np.array(results)
    print "Times, results for sim ",sim.Name()
    print times, results
    return times, results

def run(func, name, ylabel,ylog=True,compares={},ylim=None):
    global sims
    Myr = 3.15569e13
    plt.clf()
    for simname, colline in sims.iteritems():
        # Plot simulation results
        sim = Hamu.Simulation(simname)
        times, results = funcovertime(sim,func)
        # DISABLED HERE TO ALLOW ALL SIMULATIONS TO LOOP AND CACHE
        #plt.plot(times, results,colline,label=simname)
    timetobreakreleasesdoveintodawnsky()
    for simname, colline in sims.iteritems():
        # Plot lines to compare to
        ncol = 1
        sim = Hamu.Simulation(simname)
        for clabel, comp in compares.iteritems():
            ncol += 1
            newcol = colconv.to_rgb(colline[0]) # THIS SHOULD BE A COLOUR
            r,g,b = newcol
            newcol = (min(r+0.25*(ncol-1),1.0),
                      min(g+0.25*(ncol-1),1.0),
                      min(b+0.25*(ncol-1),1.0))
            newline = colline[1:] # copy
            newline = newline
            if len(newline) == 0:
                newline = "-"
            ctimes, cline = comp(sim)
            plt.plot(ctimes,cline,color=newcol,
                    linestyle=newline,
                    label=clabel+" "+simname)
    plt.xlabel(r"Time / Myr")
    plt.ylabel(ylabel)
    #plt.ylabel(r"Momentum / g cm/s")
    if ylog:
        plt.yscale("log")
    if ylim:
        plt.ylim(ylim)
    plt.legend(ncol=ncol,loc="best",fontsize="xx-small")
    plt.savefig("../plots/"+name)

if __name__=="__main__":
    Hamu.Workspace("MCRT")
    Msun = "M$_{\odot}$"
    global sims
    # Set up sims + labellings
    #timeplot.sims["MCMCtestTable"] = "b-"
    #timeplot.sims["MCMCtestConst"] = "r--"
    #timeplot.sims["MCMCtestTLowC"] = "g--"
    #timeplot.sims["MCMCtestTLDly"] = "c--"
    sims["N47_M4_B02"] = "b"
    sims["N48_M4_B02"] = "r"
    sims["N49_M4_B02"] = "g"
    sims["N00_M4_B02"] = "k"
    #timeplot.sims["N48_M4_B00"] = "r--"
    # Run for each function
    run(radialstats, "radiusstats.pdf", ylabel="Ionisation Front Radius / pc",
        compares={"M02":matzner2002.Findrii},ylim=[0.1,60])
