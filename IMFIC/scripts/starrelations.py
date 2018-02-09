'''
Correlate the structure of the cloud to the SFE
Sam Geen, January 2018
'''

from startup import *

import pymses
from pymses.filters import CellsToPoints 
from pymses.utils import constants as C

import tsfe, regression

toplot = "alltimemax"
toplot = "firstmass"
# Of course every sim formed first star at first time in the IMF runs
toplot = "firsttime"
toplot = "nphotons"
toplot = "nphotonstot"
toplot = "compactness"

def findclustercompactness(snap):
    sink = sinks.FindSinks(snap)
    x = sink.x
    y = sink.y
    z = sink.z
    m = sink.mass
    mtot = np.sum(m)
    com = np.array([np.sum(x*m)/mtot,
                    np.sum(y*m)/mtot,
                    np.sum(z*m)/mtot])
    var = np.sum((x-com[0])**2+(y-com[1])**2+(z-com[2])**2) / float(len(m))
    stddev = np.sqrt(var)
    return c

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

def _findNphotons(snap,dtins=1.0):
    # Note: dt set to 1 second by default
    stellar = stellars.FindStellar(snap)
    nphotons = np.zeros(5)
    time = stellar.time
    for tcreated, mass in zip(stellar.tcreated, stellar.mass):
        age = time - tcreated
        ageins = age * Myrins
        #dtins = dt * Myrins
        nphotons += singlestar.star_radiation(mass,ageins,dtins)
    return nphotons
findNphotons = Hamu.Algorithm(_findNphotons)

def runforsim(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    # Find last TSFE
    mmax = 0.0
    tcreatedlist = {}
    if toplot == "compactness":
        snap = sim.FindAtTime(tffcloud_code)
        c = findclustercompactness(snap)
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
        if toplot == "nphotons":
            nphot = findNphotons(snap)
            nUV = np.sum(nphot[3:5])
            mmax = max(mmax,nUV)
        if toplot == "nphotonstot":
            stellar = stellars.FindStellar(snap)
            for mass, tcreated in zip(stellar.mass, stellar.tcreated):
                if not mass in tcreatedlist:
                    tcreatedlist[mass] = tcreated
    if toplot == "nphotonstot":
        lasttimeins = 7.0 * Myrins # Hard-coded to last time we have snaps for all sims
        for mass, tcreated in tcreatedlist.iteritems():
            ageins = lasttimeins - tcreated*Myrins
            nphotons = singlestar.star_radiation(mass,0.0,ageins)
            mmax += np.sum(nphotons[3:5])
    sfe = tsfe.tsfeinsnap(sim.Snapshots()[-1])
    col = linestyles.Colour(simname)
    return sfe,mmax,col

def run(simnames,plotname,toplotin=None):
    global toplot
    if toplotin is not None:
        toplot = toplotin
    plt.clf()
    tsfes = []
    Fs = []
    cols = []
    for simname in simnames:
        sfe, F, col = runforsim(simname)
        tsfes.append(sfe)
        Fs.append(F)
        cols.append(col)
        label = linestyles.Label(simname)
        plt.scatter(F,sfe*100,s=80,
                    marker="o",c=col,edgecolors="k",
                    label=label,zorder=2)
        plt.gca().annotate(r".\ "+label, (F,sfe*100),fontsize="x-small",zorder=1)
    tsfes = np.array(tsfes)
    Fs = np.array(Fs)
    #plt.scatter(Fs,tsfes*100,s=80,marker="o",c=cols,edgecolors='k')
    #xlim = plt.gca().get_xlim()
    #plt.xlim([xlim[0],xlim[1]*2.0])
    #leg = plt.legend(fontsize="x-small")
    #leg.get_frame().set_facecolor('none')
    #linestyles.Legend("imf")
    if toplot == "alltimemax":
        plt.xlabel("Most massive star / "+Msolar)
    if toplot == "firstmass":
        plt.xlabel("Mass of first star / "+Msolar)
    if toplot == "firsttime":
        plt.xlabel("Time first star formed / Myr")
    if toplot == "nphotons":
        plt.xlabel("Peak cluster UV photon emission rate / Hz")
    if toplot == "nphotonstot":
        plt.xlabel("Total number of UV photons emitted before 7 Myr")
    allstartxt = "_"+toplot
    plt.ylabel("$\%$ TSFE (final)")
    plt.xscale("log")
    plt.yscale("log")
    plotname = "../plots/starrelations"+allstartxt+".pdf"
    regression.writecoeff(np.log10(Fs),np.log10(tsfes*100),plotname.replace(".pdf","_rxy.txt"))
    plt.savefig(plotname)

if __name__=="__main__":
    run(imfsims,"imf")

