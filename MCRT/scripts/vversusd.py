'''
Plot v versus density for each cell in the sim
Sam Geen, June 2015
'''

import customplot
import matplotlib.pyplot as plt

import os
import hydrofuncs

from pymses.filters import CellsToPoints

import numpy as np

import Hamu

import prettyplotlib as ppl
from prettyplotlib import plt
from prettyplotlib import brewer2mpl
import string
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap

red_purple = brewer2mpl.get_map('RdPu', 'Sequential', 9,reverse=True).mpl_colormap

def nvhist(snap,nbins=40):
    amr = snap.amr_source(["rho","P","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    rhos = hydrofuncs.scale_by_units(snap,"rho")
    vrads = hydrofuncs.scale_by_units(snap,"vrad")
    rhos = rhos(cells)
    vrads = vrads(cells)
    rlims = [rhos.min(),rhos.max()]
    vlims = [vrads.min(),vrads.max()]
    return np.histogram2d(np.log10(rhos),vrads,bins=nbins),rlims,vlims 

nvhistHamu = Hamu.Algorithm(nvhist)

def pramhist(snap,nbins=40):
    mH = 1.66e-24
    X = 0.76
    amr = snap.amr_source(["rho","P","vel"])
    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    rhos = hydrofuncs.scale_by_units(snap,"rho")
    vrads = hydrofuncs.scale_by_units(snap,"vrad")
    rhos = rhos(cells)
    vrads = vrads(cells)*1e5 # convert to cm/s from km/s
    pram = (vrads**2 * rhos * mH/X)[vrads < 0.0]
    hist = np.histogram(np.log10(pram),bins=nbins)
    return hist
    
pramhistHamu = Hamu.Algorithm(pramhist)

def plotsnap(snap,simname):
    histtuple,rlims,vlims = nvhistHamu(snap)
    hist, xarr,yarr = histtuple
    #import pdb; pdb.set_trace()
    xl,yl = hist.shape
    #xarr = np.arange(-rlims[0],rlims[1]*1.0000001,(rlims[1]-rlims[0])/(xl-1.0))
    #yarr = np.arange(-vlims[0],vlims[1]*1.0000001,(vlims[1]-vlims[0])/(yl-1.0))
    plt.pcolormesh(xarr,yarr,np.log10(hist.T+1),cmap=red_purple)
    plt.plot(xarr,yarr*0.0,"m:")
    #plt.xlim(rlims)
    #plt.ylim(vlims)
    plt.xlim(xarr.min(),xarr.max())
    plt.ylim(yarr.min(),yarr.max())
    try:
        os.mkdir("../plots/vversusd/"+simname)
    except:
        pass
    outnum = str(snap.OutputNumber()).zfill(5)
    plt.xlabel("log($n_H$ / cm$^{-3}$)")
    plt.ylabel("Radial Velocity / km/s")
    #plt.xscale("log")
    plt.savefig("../plots/vversusd/"+simname+"/vversusd_"+outnum+".png")
    plt.clf()
    hist,bins = pramhistHamu(snap)
    hist = np.cumsum(hist[::-1])[::-1] / float(np.sum(hist))
    #import pdb; pdb.set_trace()
    plt.plot(0.5*(bins[1:]+bins[:-1]),hist)
    plt.yscale("log")
    plt.xlabel("log(P$_{ram}$ / ergs / cm^3)")
    plt.ylabel("P( > P$_{ram}$)")
    plt.savefig("../plots/vversusd/"+simname+"/pramhist_"+outnum+".pdf")


if __name__=="__main__":
    sim = Hamu.Simulation("N00_M4_B02_C")
    snap = sim.Snapshots()[10]
    plotsnap(snap,sim.Name())
