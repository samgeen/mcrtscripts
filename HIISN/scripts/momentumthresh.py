'''
Plot properties in simulations
Sam Geen, November 2015
'''

import os, glob
import Hamu
import pymses
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import linestyles
import findproperties

gamma = 1.4 # Ramses default
mH = 1.66e-24
kB = 1.38e-16

tsn = 5.52990445



def momentuminsim(sim,nHthresh=0.0,timerange=None):
    times = []
    hfunc = Hamu.Algorithm(findproperties.momentuminsnap)
    Ts = [] # Allows generic formats for non-scalar quantities, etc
    if timerange is None:
        timerange = [-1e30,1e30]
    for snap in sim.Snapshots():
        psnap = snap.RawData()
        time = psnap.info["time"]*psnap.info["unit_time"].express(C.Myr)
        if time >= timerange[0] and time <= timerange[1]:
            if nHthresh > 0.0:
                T = hfunc(snap,nHthresh)
            else:
                T = hfunc(snap)
            Ts.append(T)
            times.append(time)
    times = np.array(times)
    # Don't cast Ts in case it's a tuple or whatever
    return times, Ts

def run(flux,linelabels=[],timerange=None):
    plt.clf()
    
    sns = ["NSN","SN"]
    simnames = ["N"+flux+"-"+x for x in sns]

    cols = ["#fd8d3c","#f03b20","#bd0026","#000000"][::-1]
    denses = [0.0,10.0,100.0,1000.0]

    plotlines = []
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        line = linestyles.line(simname)
        print line
        for col, dens in zip(cols, denses):
            times, Ts = momentuminsim(sim,dens,timerange)
            # HACK - SUBTRACT tsn
            times -= tsn
            plt.plot(times,Ts,
                     linestyle=line,
                     color=col)

        plotlines.append(mlines.Line2D([], [], 
                                       color="k",
                                       linestyle=line,
                                       label=sim.Label()))
    labels = [line.get_label() for line in plotlines]
    leg0 = plt.legend(plotlines,labels,loc="upper right",
                      ncol=1,fontsize="small",frameon=False)

    # Density limits
    #cols = ["#fecc5c","#fd8d3c","#f03b20","#bd0026"]
    plotlines = []
    denses = ["Total",
              "$n > 10 $cm$^{-3}$",
              "$n > 100 $cm$^{-3}$",
              "$n > 1000 $cm$^{-3}$"]
    for col, dens in zip(cols,denses):
        plotlines.append(mlines.Line2D([], [], 
                                       color=col, linestyle="-",
                                       label=dens))
    labels = [line.get_label() for line in plotlines]
    leg1 = plt.legend(plotlines,labels,ncol=1,fontsize="small",
                      frameon=False,loc="upper left")
    plt.gca().add_artist(leg0)
    # Quantity-specific stuff
    #plt.xscale("log")
    plt.xlabel("Time After SN / Myr")
    plt.yscale("log")
    plt.ylabel("momentum / g cm/s")
    #plt.ylim([4e43,2e45])
    plt.savefig("../plots/momentumthresh_N"+flux+".pdf")

if __name__=="__main__":
    fluxes = ["00","49","50","51"]
    tr = [tsn-0.1,tsn+3.0]
    for flux in fluxes:
        # HACK TO FORCE THE SIMULATIONS TO RUN
        run(flux,timerange=tr)

