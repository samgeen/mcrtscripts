'''
Correlate the structure of the cloud to the SFE
Sam Geen, January 2018
'''

from startup import *

import pymses
from pymses.filters import CellsToPoints 
from pymses.utils import constants as C

import tsfe, regression, trackstars, findproperties

# List of possible plotting values
#toplot = "alltimemax"
#toplot = "firstmass"
# Of course every sim formed first star at first time in the IMF runs
#toplot = "firsttime"
#toplot = "nphotons"
#toplot = "nphotonstot"
#toplot = "compactness"
#toplot = "tracklength"
#toplot = "nphotonstff"

def findclustercompactness(snap):
    sink = sinks.FindSinks(snap)
    x = sink.x
    y = sink.y
    z = sink.z
    m = sink.mass
    if len(m) <= 1:
        return None
    mtot = np.sum(m)
    com = np.array([np.sum(x*m)/mtot,
                    np.sum(y*m)/mtot,
                    np.sum(z*m)/mtot])
    var = np.sum((x-com[0])**2+(y-com[1])**2+(z-com[2])**2) / float(len(m))
    stddev = np.sqrt(var)
    return stddev

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

def runforsim(simname,xvalue,yvalue="sfe"):
    print "Running for simulation", simname
    sim = hamusims[simname]
    FindStellarHamu = Hamu.Algorithm(stellars.FindStellar)
    xout = 0.0
    tcreatedlist = {}
    toplot = xvalue
    if toplot == "compactness":
        snap = sim.FindAtTime(1.5*tffcloud_code)
        c = findclustercompactness(snap)
        xout = c
    if toplot == "tracklength":
        tracks = trackstars.runforsim(simname)
        tracklength = 0.0
        for track in tracks.itervalues():
            tracklength += np.sqrt(np.sum((track[0,:]-track[track.shape[0]-1,:])**2))
        tracklength /= float(len(tracks.keys()))
        xout = tracklength
    for snap in sim.Snapshots():
        if toplot == "alltimemax":
            xout = max(findmaxstar(snap),xout)
        if toplot == "firstmass":
            xout = findmaxstar(snap)
            if xout > 0.0:
                break
        if toplot == "firsttime":
            tcreated = findtcreated(snap)
            if tcreated > 0.0:
                xout = tcreated
                break
        if toplot == "nphotons":
            nphot = findNphotons(snap)
            nUV = np.sum(nphot[3:5])
            xout = max(xout,nUV)
        if toplot == "nphotonstot" or toplot == "nphotonstff":
            stellar = FindStellarHamu(snap)
            for mass, tcreated in zip(stellar.mass, stellar.tcreated):
                if not mass in tcreatedlist:
                    tcreatedlist[mass] = tcreated
    if toplot == "nphotonstot":
        lasttimeins = 7.0 * Myrins # Hard-coded to last time we have snaps for all sims
        for mass, tcreated in tcreatedlist.iteritems():
            ageins = lasttimeins - tcreated*Myrins
            nphotons = singlestar.star_radiation(mass,0.0,ageins)
            xout += np.sum(nphotons[3:5])
    if toplot == "nphotonstff":
        lasttimeins = tffcloud_Myr * Myrins
        for mass, tcreated in tcreatedlist.iteritems():
            ageins = lasttimeins - tcreated*Myrins
            nphotons = singlestar.star_radiation(mass,0.0,ageins)
            xout += np.sum(nphotons[3:5])
    if yvalue == "sfe":
        yout = tsfe.tsfeinsnap(sim.Snapshots()[-1]) * 100 # as a percentage
    if yvalue == "momentum":
        func = Hamu.Algorithm(findproperties.totalmomentuminsnap)
        yout = func(sim.FindAtTime(tffcloud_code))
    if xvalue == "momentum":
        func = Hamu.Algorithm(findproperties.totalmomentuminsnap)
        xout = func(sim.FindAtTime(tffcloud_code))
    rfunc = Hamu.Algorithm(findproperties.radiusinsnap)
    if xvalue == "radius":
        xout = rfunc(sim.FindAtTime(tffcloud_code))
    if yvalue == "radius":
        yout = rfunc(sim.FindAtTime(tffcloud_code))
    return xout,yout

def run(simnames,plotname,xvalue,yvalue="sfe"):
    plt.clf()
    xvals = []
    yvals = []
    for simname in simnames:
        xval, yval = runforsim(simname,xvalue,yvalue)
        col = linestyles.Colour(simname)
        if xval is not None:
            xvals.append(xval)
            yvals.append(yval)
            label = linestyles.Label(simname)
            plt.scatter(xval,yval,s=80,
                        marker="o",c=col,edgecolors="k",
                        label=label,zorder=2)
            plt.gca().annotate(". "+label, (xval,yval),fontsize="x-small",zorder=1,clip_on=True)
    xvals = np.array(xvals)
    yvals = np.array(yvals)
    # Set xlabel
    toplot = xvalue
    if toplot == "alltimemax":
        plt.xlabel("Most massive star / "+Msolar)
    if toplot == "firstmass":
        plt.xlabel("Mass of first star / "+Msolar)
    if toplot == "firsttime":
        plt.xlabel("Time first star formed / Myr")
    if toplot == "nphotons":
        plt.xlabel("Peak cluster UV photon emission rate / Hz")
    if toplot == "nphotonstot":
        plt.xlabel("Total number of UV photons before 7 Myr")
    if toplot == "nphotonstff":
        plt.xlabel(r"Total number of UV photons in 1 freefall time")
    if toplot == "compactness":
        plt.xlabel("Stddev of sink separation / pc")
    if toplot == "tracklength":
        plt.xlabel("Average distance a massive star travels / pc")
    allstartxt = "_"+toplot
    plt.xscale("log")
    plt.yscale("log")
    momtxt = "Momentum / g cm / s"
    if xvalue == "momentum":
        plt.xlabel(momtxt)
    if yvalue == "momentum":
        plt.ylabel(momtxt)
    if yvalue == "sfe":
        plt.ylabel("$\%$ TSFE (final)")
    rtxt = "Mean radius of HII region / pc"
    if yvalue == "radius":
        plt.ylabel(rtxt)
    if xvalue == "radius":
        plt.xlabel(rtxt)
    if yvalue != "sfe":
        allstartxt += "_"+yvalue
    # Stop outflowing labels ruining the plot area
    xlim = plt.gca().get_xlim()
    plt.xlim([xlim[0],xlim[1]*1.9])
    plt.ylim([5.0,22.0])
    plotname = "../plots/starrelations"+allstartxt+plotname+".pdf"
    regression.writecoeff(np.log10(xvals),np.log10(yvals),plotname.replace(".pdf","_rxy.txt"))
    #plt.gcf().set_size_inches(18.5, 10.5)
    plt.savefig(plotname,bbox_inches='tight',dpi=150)

if __name__=="__main__":
    run(imfsims,"imf","nphotonstff","sfe")

