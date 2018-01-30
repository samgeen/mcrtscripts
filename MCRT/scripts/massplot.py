'''
Plots the various masses in the simulation
Sam Geen, March 2015
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import Hamu

import timeplot

import pymses
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

def FindTotal(cells, snap):
    masses = cells["rho"]*cells.get_sizes()**3.0
    masses *= snap.info["unit_mass"].express(C.Msun)
    return np.nansum(masses) # Total mass

def FindIonised(cells, snap):
    masses = cells["rho"]*cells.get_sizes()**3.0
    masses *= snap.info["unit_mass"].express(C.Msun)
    masses *= cells["xHII"] # mass in ions = mass in cell * ionisation fraction
    return np.nansum(masses)

def FindNeutral(cells, snap):
    masses = cells["rho"]*cells.get_sizes()**3.0
    masses *= snap.info["unit_mass"].express(C.Msun)
    masses *= (1.0-cells["xHII"])
    return np.nansum(masses)

def FindUnstable(cells, snap):
    '''
    Find cells that are Jeans unstable
    '''
    kB = 1.3806488e-16
    rhos = cells["rho"]*snap.info["unit_density"].express(C.H_cc)
    cellsizes = cells.get_sizes()*snap.info["unit_length"].express(C.pc)
    T = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    cs = np.sqrt(1.4*kB*T/1.6/1.66e-24)/1e5 # in km/s
    nsf = 1e3 * ((cellsizes/0.4) / (cs/0.2))**(-2.0)
    masses = cells["rho"]*cells.get_sizes()**3.0
    masses *= snap.info["unit_mass"].express(C.Msun)
    #rlim = 2912.0 # cm-3; given cs = 0.1km/s, l_jeans = 120pc/2^10
    sfcells = rhos > nsf
    try:
        mass = np.sum(masses[sfcells])
    except:
        mass = 0.0
    return mass

def MassesInSnap(snap):
    '''
    Find the various masses in the snapshot
    '''
    amr = snap.amr_source(["rho","P","xHII"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    unstable = FindUnstable(cells, snap)
    ionised  = FindIonised(cells, snap)
    total    = FindTotal(cells, snap)
    neutral  = total - ionised # FindNeutral(cells, snap)
    return (unstable, ionised, neutral, total)

def plotforsim(times, results,col,label):
    '''
    Plot for a given simulation
    results is array of shape len(times) x 4
    '''
    # Unstable
    plt.plot(times,results[:,0],col+":")
    # Neutral
    plt.plot(times,results[:,2],col+"-.")
    # Ionised
    plt.plot(times,results[:,1],col+"--")
    # Total
    plt.plot(times,results[:,3],col+"-",label=label)

def run():
    sims = timeplot.sims
    folder = timeplot.folder
    if folder:
        try:
            os.mkdir("../plots/"+folder)
        except:
            pass # Eh
    Myr = 3.15569e13
    plt.clf()
    for simname, colline in sims.iteritems():
        # Plot simulation results
        sim = Hamu.Simulation(simname)
        times, results = timeplot.funcovertime(sim,MassesInSnap)
        totalmass = results[0,3]
        results /= totalmass # Do mass fraction
        plotforsim(times, results, colline[0],sim.Label())
    '''
    firstc = True
    for simname, colline in sims.iteritems():
        # Plot lines to compare to
        ncol = 1
        sim = Hamu.Simulation(simname)
        tstart = 1.25
        try:
            tstart = starts[sim.Name()]
        except:
            print "Didn't find tstart in timeplot.run, using default 1.25Myr"
            print sim.Name(), starts.keys()
        matzner2002.tstart = tstart
        hennebellemodel.tstart = tstart
        outflowmodel.tstart = tstart
        lines = ["-","-","--",":"]
        if "momentum" in name:
            outflowmodel.momback = FindMomBack(sim, starts[simname])
        for clabel, comp in compares.iteritems():
            ncol += 1
            #newcol = colconv.to_rgb(colline[0]) # THIS SHOULD BE A COLOUR
            #r,g,b = newcol
            #newcol = (min(r+0.333*(ncol-1),1.0),
            #          min(g+0.333*(ncol-1),1.0),
            #          min(b+0.333*(ncol-1),1.0))
            #newline = colline[1:] # copy
            #newline = newline
            # HACK - now use dots and dashes for compares
            newcol = colline[0]
            newline = lines[ncol]
            #if len(newline) == 0:
            #    newline = "-"
            ctimes, cline = comp(sim)
            currlabel = clabel
            if not firstc:
                currlabel = ""
            plt.plot(ctimes,cline,color=newcol,
                    linestyle=newline,
                    label=currlabel)
        firstc = False
    '''
    plt.xlabel(r"Time / Myr")
    #plt.ylabel("Mass / M$_{\odot}$")
    plt.ylabel(r"Mass Fraction of Cloud")
    plt.xlim([0,5])
    plt.ylim([1e-3,1.0])
    #plt.ylabel(r"Momentum / g cm/s")
    plt.yscale("log")
    simleg = plt.legend(loc="best",fontsize="xx-small")
    # Mass type legen
    uline     = mlines.Line2D([],[],linestyle=":",color="k")
    nline    = mlines.Line2D([],[],linestyle="-.",color="k")
    iline    = mlines.Line2D([],[],linestyle="--",color="k")
    tline   = mlines.Line2D([],[],linestyle="-",color="k")
    massleg = plt.legend((uline,nline,iline,tline),
                         ("Unstable", "Neutral", "Ionised",
                          "Total"),fontsize="xx-small")
    plt.gca().add_artist(massleg)
    plt.gca().add_artist(simleg)
    plt.savefig("../plots/"+folder+"/masses.pdf")
