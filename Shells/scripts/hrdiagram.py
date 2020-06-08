"""
Make a HR diagram out of the stellar evolution tables
Sam Geen, May 2020
"""

import matplotlib.pyplot as plt

import numpy as np

import singlestar, readtefftables

_tablesetup = False
_singlestarLocation = ""

def moreatstartandend(x):
    # Gives more values at start and end
    # x must be 0 to 1
    return 0.5 - 0.5*np.cos(np.pi*x)

def readtable(metal,rotating):
    global _tablesetup
    global _singlestarlocation
    rtext = ""
    if not rotating:
        rtext = "_notrotating"
    singlestarLocation = "/home/samgeen/Programming/Astro/StellarSources/Outputs/singlestar_z"+str(metal)+rtext
    if _tablesetup and singlestarLocation != _singlestarLocation:
        singlestar.star_reset()
    singlestar.star_setup(singlestarLocation)
    _tablesetup = True

def plotformetal(masses,metal,rotating,iband):
    print "Plotting for metal", metal, "rotating:", rotating
    readtable(metal,rotating)
    print "Read table"
    times = np.linspace(0,1.0,10000)
    if rotating:
        line = "-"
    else:
        line = "--"
    if metal == 0.014:
        cmap = plt.get_cmap("Reds")
    else:
        cmap = plt.get_cmap("Blues")
    for imass, mass in enumerate(masses):
        colour = cmap(0.5+0.5*float(imass)/float(len(masses)))
        lifetime = singlestar.star_lifetime(mass)
        startimes = moreatstartandend(times)*lifetime
        # Kband, Vband, Lbol 
        lums = []
        Teffs = readtefftables.Teff(startimes,lifetime,mass,metal,rotating).flatten()
        for time in startimes:
            dt = 1e7
            bands = singlestar.star_bandenergies(mass,time,dt)/dt
            lums.append(bands[iband])
        lums = np.array(lums)
        mask = lums*Teffs > 1e20
        lums = lums[mask]
        Teffs = Teffs[mask]
        plt.plot(Teffs,lums,linestyle=line,color=colour)
    print "Done!"


def plot(band):
    print "---"
    print "PLOTTING FOR BAND", band
    plt.clf()
    if band == "LK":
        ylabel = "$L_{K}$"
        iband = 0
    if band == "LV":
        ylabel = "$L_{V}$"
        iband = 1
    if band == "Lbol":
        ylabel = "$L_{bol}$"
        iband = 2
    masses = [15,30,60,120]
    # Loop through values
    for metal in [0.014]:
        for rotating in [True,False]:
            plotformetal(masses,metal,rotating,iband)
    # Do useful plot stuff
    plt.xlabel("Temperature / K")
    plt.ylabel(ylabel + " / erg/s")
    #plt.xscale("log")
    plt.yscale("log")
    # Do weird inverted x axis thing for HR diagrams
    plt.gca().invert_xaxis()
    # Save
    plt.savefig("../plots/hrdiagram_Tvs"+band+".pdf")

if __name__=="__main__":
    [plot(band) for band in ["LK","LV","Lbol"]]

