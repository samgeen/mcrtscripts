'''
Model the star-forming mass in the cloud very briefly
Sam Geen, March 2015
'''

import customplot
import matplotlib.pyplot as plt

import matplotlib.lines as mlines

import Hamu

import numpy as np

import os, errno

import pymses
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import findclumps
import massplot, timeplot
# This is necessary for cache file loading of ClumpMaker objects
from findclumps import ClumpMaker

import plottimeseries

import linestyles

def Accretion(minit,times,tff):
    # TODO: Try exponential growth
    rate = minit / tff
    masses = minit + minit* times / tff # np.exp(times / tff)
    return masses

def PlotAccretion(simnames,name,nomodel=False):
    plt.clf()
    #cols = ["k","b","m","r"]
    cols = []
    for simname in simnames:
        cols.append(linestyles.col(simname))
    icol = 0
    tstart = 1.25
    tff_dense = 1.45
    myr = None
    tff_diffuse = 4.4 # Myr - calculated at 1e2 at/cc, mu = 1.4
    sims = {}
    tstarts = {}
    tff_denses = {}
    tff_diffuses = {}
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        sims[simname] = sim
        if "_C2" in simname:
            tstarts[simname] = tstart * 0.75**3
            tff_denses[simname] = tff_dense * 0.75**3
            tff_diffuses[simname] = 2.0 # Myr - calculated at 5e2 at/cc
        elif "_C" in simname:
            tstarts[simname] = tstart * 0.5**3
            tff_denses[simname] = tff_dense * 0.5**3
            tff_diffuses[simname] = 43.7 # Myr Only the 1 atcc medium outside
        else:
            tstarts[simname] = tstart
            tff_denses[simname] = tff_dense
            tff_diffuses[simname] = tff_diffuse
    #sim = Hamu.Simulation("N00_M4_B02")
    sffunc = Hamu.Algorithm(plottimeseries.sfdensityinsnap)
    icol = 0
    unstable_nophotons = None
    photsnames = ["None","10$^{47}$","10$^{48}$","10$^{49}$"]
    simstomodel = simnames
    if nomodel:
        simstomodel = []
    for simname in simnames:
        sim = sims[simname]
        col = cols[icol]
        times = []
        sfmasses = []
        unionmass = []
        for snap in sim.Snapshots():
            if myr is None:
                myr = snap.RawData().info["unit_time"].express(C.Myr)
            #sfmass = sffunc(snap)
            times.append(snap.Time()*myr)
            
            #unstable, ion, neut, tot = massplot.MassesInSnapHamu(snap)
            #sfmasses.append(unstable) # unstable mass
            #unionmass.append(tot-ion) # unionised mass
        # Get masses from simulation
        timeplot.starts[simname] = tstarts[simname]
        times, results = timeplot.funcovertime(sim,massplot.MassesInSnap)
        times += tstarts[simname]
        unstable = results[:,0]
        ionised = results[:,1]
        total = results[:,3]
        # HACK - 5 is the end of the accretion phase
        if "N00" in simname:
            unstable_nophotons = unstable
            #unionised = unstable * 0.0
        else:
            #unionised = unstable_nophotons[:len(ionised)] - ionised
            pass
        # Plot
        #import pdb; pdb.set_trace()
        line = "-"
        if simname == "N00_M4_B02_C" or simname == "N00_M4_B02_C2":
            line = "--"
        plt.plot(times,unstable,color=col,linestyle=line,
                 label=sim.Label())
        #plt.plot(times,unionised,col+"--",label=sim.Label()+"T-I")
        icol += 1
    icol = 0
    for simname in simstomodel:
        col = cols[icol]
        sim = sims[simname]
        sim00name = "N00_M4_B02"
        if "_C2" in simname:
            sim00name = "N00_M4_B02_C2"
        elif "_C" in simname:
            sim00name = "N00_M4_B02_C"
        try:
            sim00 = sims[sim00name]
        except:
            sim00 = Hamu.Simulation(sim00name)
        #times, masses = findclumps.EvaporationRate(simname,time=tstart,
        #                                           mag=False)
        if myr is None:
            snap = sim.Snapshots()[0]
            myr = snap.RawData().info["unit_time"].express(C.Myr)
        snap = sim.FindAtTime(tstarts[simname]/myr)
        #masses = Accretion(masses[0],times,tff_diffuse)
        # Run accretion / evaporation model
        tstart = tstarts[simname]
        tff_dense = tff_denses[simname]
        tff_diffuse = tff_diffuses[simname]
        times, masses = findclumps.AccretionEvaporation(sim,tff_diffuse,
                                                        tstart=tff_dense,
                                                        mag=False,
                                                        nofluxsim=sim00)
        times += tff_dense#tstart
        plt.plot(times, masses,color=col,linestyle=":")
                 #label="Model")
        icol += 1

    ncol = 2
    #if nomodel:
    #    ncol = 1
    if name == "compact":
        ncol = 1
    leg1 = plt.legend(fontsize="small",loc="lower right",
                     ncol=ncol, frameon=False)

    sline = mlines.Line2D([], [], color='k', label='Simulation')
    mline = mlines.Line2D([], [], color='k', linestyle=":",
                              label='Model')
    if name == "photons":
        leg2 = plt.legend(handles=[sline,mline],ncol=2,fontsize="small",
                          frameon=False,loc="upper left")
        plt.gca().add_artist(leg1)

    # Re-align the y-axes THANK YOU VERY MUCH FOR DOING THIS YOURSELF PYPLOT
    #if not nomodel:
    #    txts = leg.get_texts()
    #    nsim = len(sims)
    #    for itxt in range(0,nsim):
    #        left = txts[itxt]
    #        right = txts[itxt+nsim]
    #        right.set_y(left.get_position()[1])
    #        left.set_verticalalignment("baseline")
    #        right.set_verticalalignment("baseline")

    plt.yscale("log")
    plt.ylabel("Total Jeans-Unstable Mass / M$_{\odot}$")
    plt.xlabel("Time / Myr")
    plt.ylim([1e3,1e4])
    nomodeltext = ""
    if nomodel:
        nomodeltext = "_nomodel"
    plt.savefig("../plots/clumps/"+name+"/sfmass"+nomodeltext+".pdf")
    print simnames

def Plot(simnames):
    plt.clf()
    icol = 0
    for simname in simnames:
        col = cols[icol]
        sim = Hamu.Simulation(simname)
        times, masses = findclumps.EvaporationRate(simname,time=1.25,mag=False)
        plt.plot(times, masses,col,
                 label=sim.Label()+" Thermal")
        icol += 1
    plt.legend(fontsize="x-small",loc="upper right")
    plt.yscale("log")
    plt.ylabel("Total Mass / M$_{\odot}$")
    plt.xlabel("Time / Myr")
    plt.savefig("../plots/clumps/photons/evaporation.pdf")
    print simnames

if __name__=="__main__":
    simnames = ["N00_M4_B02","N47_M4_B02","N48_M4_B02",
                "N48_M4_B00","N49_M4_B02"]
    #PlotAccretion(simnames,"photons",nomodel=True)
    PlotAccretion(simnames,"photons")
    #simnames = ["N00_M4_B02","N48_M4_B02",
    #            "N00_M4_B02_C2","N48_M4_B02_C2",
    #            "N00_M4_B02_C","N48_M4_B02_C"]
    simnames = ["N00_M4_B02","N48_M4_B02",
                "N00_M4_B02_C","N48_M4_B02_C",
                "N00_M4_B02_C2","N48_M4_B02_C2"]
    PlotAccretion(simnames,"compact")

