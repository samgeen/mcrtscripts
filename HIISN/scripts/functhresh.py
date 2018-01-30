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

tsn = 0.0

class MomentumFunc(object):
    '''
    Functionoid to calculate momentum
    '''

    name = "momentum"
    func = [findproperties.momentuminsnap]
    def __init__(self):
        pass

    def setupplot(self):
        plt.xlabel("Time After SN / Myr")
        plt.ylim([7e41,4e44])
        plt.yscale("log")
        plt.ylabel("momentum / g cm/s")

class MassFunc(object):
    '''
    Functionoid to calculate mass above a given threshold
    '''

    name = "mass"
    func = [findproperties.massinsnap]
    def __init__(self):
        pass

    def setupplot(self):
        plt.xlabel("Time After SN / Myr")
        plt.yscale("log")
        plt.ylabel("Mass / M$_{\odot}$")

def funcinsim(sim,f,nHlow=0.0,nHhigh=0.0,timerange=None):
    times = []
    # [0] because the func is in a list to prevent binding it as a member
    #    function to the object func
    hfunc = Hamu.Algorithm(f.func[0])
    Ts = [] # Allows generic formats for non-scalar quantities, etc
    if timerange is None:
        timerange = [-1e30,1e30]
    for snap in sim.Snapshots():
        psnap = snap.RawData()
        time = psnap.info["time"]*psnap.info["unit_time"].express(C.Myr)
        if time >= timerange[0] and time <= timerange[1]:
            if nHlow > 0.0 or nHhigh > 0.0:
                T = hfunc(snap,nHlow,nHhigh)
            else:
                T = hfunc(snap)
            Ts.append(T)
            times.append(time)
    times = np.array(times)
    # Don't cast Ts in case it's a tuple or whatever
    return times, Ts

def run(f,flux,linelabels=[],timerange=None,lowerbounding=True):
    '''
    lowerbounding - if set, only give values inside each threshold
                  - if unset, give values for all below the threshold
    '''
    plt.clf()
    
    sns = ["NSN","SN"]
    simnames = ["N"+flux+"-"+x for x in sns]

    cols = ["#fd8d3c","#f03b20","#bd0026","#000000"][::-1]
    lows = [0.0,0.0,10.0,100.0]
    highs = [0.0,10.0,100.0,1e30]

    plotlines = []
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        line = linestyles.line(simname)
        print line
        for col, low, high in zip(cols, lows, highs):
            times, Ts = funcinsim(sim,f,low,high,timerange)
            Ts = np.array(Ts)
            toplot = Ts
            # HACK - SUBTRACT tsn
            times -= tsn
            plt.plot(times,toplot,
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
              "$n < 10 $cm$^{-3}$",
              "10 cm$^{3} < n < 100 $cm$^{-3}$",
              "100 cm$^{3} < n$"] # < 1000 $cm$^{-3}$",]
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
    f.setupplot()
    #plt.ylim([4e43,2e45])
    ws = Hamu.CurrentWorkspace()
    suffix = ""
    if lowerbounding:
        suffix += "lb"
    plt.savefig("../plots/"+f.name+"thresh_N"+flux+"_"+ws+"_"+suffix+".pdf")

if __name__=="__main__":
    workspaces = ["HIISN","HIISN4"]
    #workspaces = ["HIISN4"]
    for ws in workspaces:
        Hamu.Workspace(ws)
        if "4" in ws:
            tsn = 4.25
            fluxes = ["00","47","48","49"]
        else:
            tsn = 5.52990445
            fluxes = ["00","49","50","51"]
        tr = [tsn-0.1,tsn+3.0]
        #funcs = [MassFunc(),MomentumFunc()]
        funcs = [MomentumFunc()]
        for f in funcs:
            for flux in fluxes:
                # HACK TO FORCE THE SIMULATIONS TO RUN
                run(f,flux,timerange=tr,lowerbounding=True)

