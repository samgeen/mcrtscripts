'''
Correlate the structure of the cloud to the SFE
Sam Geen, January 2018
'''

from startup import *

import pymses
from pymses.filters import CellsToPoints 
from pymses.utils import constants as C

import tsfe

toplot = "alltimemax"
toplot = "firstmass"
# Of course every sim formed first star at first time in the IMF runs
toplot = "firsttime"

def _findmaxstar(snap):
    stellar = stellars.FindStellar(snap)
    try:
        mmax = stellar.mass.max()
    except:
        mmax = 0.0
    return mmax
findmaxstar = Hamu.Algorithm(_findmaxstar)

def _findtcreated(snap):
    stellar = stellars.FindStellar(snap)
    try:
        tcreated = stellar.tcreated.min()
    except:
        tcreated = 0.0
    return tcreated
findtcreated = Hamu.Algorithm(_findtcreated)

def runforsim(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    # Find last TSFE
    mmax = 0.0
    for snap in sim.Snapshots():
        if toplot == "alltimemax":
            mmax = max(findmaxstar(snap),mmax)
        if toplot == "firstmass":
            mmax = findmaxstar(snap)
            if mmax > 0.0:
                break
        if toplot == "firsttime":
            tcreated = findtcreated(snap)
            if tcreated > 0.0:
                mmax = tcreated
                break
    sfe = tsfe.tsfeinsnap(sim.Snapshots()[-1])
    col = linestyles.Colour(simname)
    return sfe,mmax,col

def run(simnames,plotname):
    plt.clf()
    tsfes = []
    Fs = []
    cols = []
    for simname in simnames:
        sfe, F, col = runforsim(simname)
        tsfes.append(sfe)
        Fs.append(F)
        cols.append(col)
    tsfes = np.array(tsfes)
    Fs = np.array(Fs)
    plt.scatter(Fs,tsfes*100,s=80,marker="o",c=cols,edgecolors='k')
    print tsfes
    if toplot == "alltimemax":
        plt.xlabel("Most massive star / "+Msolar)
    if toplot == "firstmass":
        plt.xlabel("Mass of first star / "+Msolar)
    if toplot == "firsttime":
        plt.xlabel("Time first star formed / Myr")
    allstartxt = "_"+toplot
    plt.ylabel("$\%$ TSFE (final)")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("../plots/starrelations"+allstartxt+".pdf")

if __name__=="__main__":
    run(imfsims,"imf")

