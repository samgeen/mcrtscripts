'''
Plot the ionised masses in each run to test convergence
Sam Geen, February 2014
'''

import os
import Hamu
import numpy as np
import matplotlib.pyplot as plt

import pymses
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

def massfractioninsnap(snap):
    '''
    Find the cumulative probability per mass
    '''
    amr = snap.amr_source(["rho"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    rhos = cells["rho"]*snap.info["unit_density"].express(C.H_cc)
    ddens = rhos
    vols = cells.get_sizes()**3.0
    mass = vols * cells["rho"]*snap.info["unit_mass"].express(C.Msun)
    inds = np.argsort(rhos)
    mass = mass[inds]
    rhos = rhos[inds]
    cmass = np.cumsum(mass)
    # Return cumulative mass as a function of density
    return np.log10(rhos), cmass
    '''
    pdf, edges = np.histogram(ddens,bins=50,weights=mass)
    pdf = pdf.astype(np.float32)
    pdf /= np.sum(mass)
    # Turn bin edges into bin centre
    dens = 0.5 * (edges[0:len(edges)-1] + edges[1:len(edges)])
    print dens, pdf #edges[0:len(edges)-1] + edges[1:len(edges)]
    return dens, pdf
    '''

def densitypdfinsnap(snap):
    #h = snap.hamusnap
    #path = h.Path()
    #path = path[0:path.find("output_")]
    #ro = pymses.RamsesOutput(path,h.OutputNumber())
    amr = snap.amr_source(["rho"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    rhos = cells["rho"]*snap.info["unit_density"].express(C.H_cc)
    ddens = rhos
    vols = cells.get_sizes()**3.0
    pdf, edges = np.histogram(np.log10(ddens),bins=50,weights=vols)
    pdf = pdf.astype(np.float32)
    pdf /= np.sum(vols)
    # Turn bin edges into bin centre
    dens = 0.5 * (edges[0:len(edges)-1] + edges[1:len(edges)])
    print edges[0:len(edges)-1] - edges[1:len(edges)]
    return dens, pdf

masspdfHamu = Hamu.Algorithm(massfractioninsnap)
densitypdfHamu = Hamu.Algorithm(densitypdfinsnap)

def makeplot(timeIn = 4.0):
    Myr = 3.15569e13
    timeToSample = timeIn+1.25 # 1.25 is the starting time
    plt.clf()
    ns = ["49","48","47","00"]
    bs = ["02"]
    cols = ["r","g","b","k"]
    lines = ["-"]
    times = list()
    for b, line in zip(bs, lines):
        for n, col in zip(ns, cols):
            simname = "N"+n+"_M4_B"+b
            sim = Hamu.Simulation(simname)
            # Get time units
            firstSnap = sim.Snapshots()[0].RawData()
            unit_time = firstSnap.info["unit_time"].express(C.Myr)
            # Get snapshot at time
            snap = sim.FindAtTime(timeToSample/unit_time)
            tMyr = snap.Time() * unit_time
            dens, pdf = masspdfHamu(snap)
            plt.plot(dens, pdf,col+line,label=sim.Label())
            times.append(tMyr-1.25)
    #plt.xlabel(r"log(n$_H$ / cm$^{-3}$)")
    #plt.ylabel(r"P(n$_H$)")
    plt.xlabel("n$_H$")
    plt.ylabel("Cumulative mass")
    plt.yscale("log")
    plt.ylim([1e1,1e5])
    #plt.xlim([-2,4])
    #plt.ylim([5e-4,1])
    plt.legend(ncol=2)
    plt.savefig("../plots/densitypdf/densitypdf_"+str(int(tMyr-1.25))+".pdf")
    print times

def run():
    for t in [2.0]: #np.arange(0.0,6.0,1.0):
        makeplot(t)

if __name__=="__main__":
    run()
