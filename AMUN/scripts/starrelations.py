'''
Correlate the structure of the cloud to the SFE
Sam Geen, January 2018
'''

from startup import *

import pymses
from pymses.filters import CellsToPoints 
from pymses.utils import constants as C

import regression, trackstars, findproperties, plotproperties

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

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
        yout = plotproperties.tsfeinsnap(sim.Snapshots()[-1]) * 100 # as a percentage
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

def makeplot(ax,simnames,plotname,xvalue,yvalue="sfe",yon=True):
    xvals = []
    yvals = []
    for simname in simnames:
        xval, yval = runforsim(simname,xvalue,yvalue)
        col = linestyles.Colour(simname)
        if xval is not None:
            xvals.append(xval)
            yvals.append(yval)
            label = linestyles.Label(simname)
            ax.scatter(xval,yval,s=80,
                        marker="o",c=col,edgecolors="k",
                        label=label,zorder=2)
            ax.annotate(". "+label, (xval,yval),fontsize="x-small",zorder=1,clip_on=True)
    xvals = np.array(xvals)
    yvals = np.array(yvals)
    # Set xlabel
    toplot = xvalue
    plotminor = True
    plotnox = False
    if toplot == "alltimemax":
        xtxt = "Most massive star $M_{max}$ / "+Msolar
        ax.set_xlim([5,120])
        ax.get_xaxis().set_major_locator(MultipleLocator(20))
        #ax.get_xaxis().set_minor_locator(MultipleLocator(20000))
        plotminor = False
        #ax.tick_params(axis='x',which='minor',bottom='off')
    if toplot == "firstmass":
        xtxt = "Mass of first star $M_{first}$ / "+Msolar
        ax.set_xlim([5,100])
        ax.get_xaxis().set_major_locator(MultipleLocator(20))
        #ax.tick_params(axis='x',which='minor',bottom='off')
        plotminor = False
    if toplot == "firsttime":
        xtxt = "Time first star formed / Myr"
    if toplot == "nphotons":
        xtxt = "Peak photon emission rate $S_{*,max}$ / Hz"
    if toplot == "nphotonstot":
        xtxt = "Total photons before 7 Myr, $N_{tot}$"
    if toplot == "nphotonstff":
        xtxt = r"Total photons in first $t_{ff}$, $N_{ff}$"
    if toplot == "compactness":
        xtxt = "Cluster Compactness $R_{RMS}$ / pc"
        xticks = [5,6,7,8,9]
        xticklabels = [str(x) for x in xticks]
        plotnox = True
    if toplot == "tracklength":
        xtxt = "Average distance a massive star travels $D_{track}$ / pc"
        xticks = [3,4,5,6,7,8,9,10]
        xticklabels = [str(x) for x in xticks]
        plotnox = True
    allstartxt = "_"+toplot
    ax.set_xscale("log")
    ax.set_yscale("log")
    momtxt = "Momentum / g cm / s"
    if xvalue == "momentum":
        xtxt = momtxt
    if yvalue == "momentum":
        ytxt = momtxt
    if yvalue == "sfe":
        ytxt = "SFE (final)"
    rtxt = "Mean radius of HII region / pc"
    if yvalue == "radius":
        ytxt = rtxt
    if xvalue == "radius":
        xtxt = rtxt
    if yvalue != "sfe":
        allstartxt += "_"+yvalue
    # Stop outflowing labels ruining the plot area
    xlim = ax.get_xlim()
    ax.set_xlim([xlim[0],xlim[1]*1.2])
    #ax.set_ylim([5.0,22.0])    
    ax.set_xlabel(xtxt)
    ticks = range(5,23)
    ticklabels = {tick:" " for tick in ticks}
    ticklabelsempty = {tick:" " for tick in ticks}
    ticklabels[5] = "5\%"
    ticklabels[10] = "10\%"
    ticklabels[20] = "20\%"
    if xlim[1] < 1000:
        if plotnox:
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.xaxis.set_minor_locator(plt.NullLocator())
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
        else:
            if plotminor:
                ax.get_xaxis().set_minor_formatter(plt.ScalarFormatter())
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        #else:
        #    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_minor_locator(plt.NullLocator())
    if yon:
        ax.set_ylabel(ytxt)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabels.values())
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabelsempty.values())
    plotname = "../plots/starrelations"+allstartxt+plotname+".pdf"
    regression.writecoeff(np.log10(xvals),np.log10(yvals),plotname.replace(".pdf","_rxy.txt"))
    #plt.gcf().set_size_inches(18.5, 10.5)
    
def run(simnames,plotname,xvalues,yvalue="sfe",suffix="all"):
    nplots = len(xvalues)
    ncols  = 3
    nrows  = nplots/ncols
    if nplots < ncols:
        ncols = nplots
        nrows = 1
    if nrows*ncols < nplots:
        nrows += 1
    fig, axes = plt.subplots(nrows,ncols)#,sharey=True)
    try:
        dum = axes.shape
    except:
        axes = np.array([axes])
    irow = 0
    icol = 0
    yon = True
    for ax, xvalue in zip(axes.flatten(), xvalues):
        makeplot(ax,simnames,plotname,xvalue,yvalue="sfe",yon=yon)
        yon = False
        icol += 1
        #if not yon:
        #    ax.get_yaxis().set_ticklabels([])
        if icol == ncols:
            icol = 0
            irow += 1
            yon = True
    fig.subplots_adjust(wspace=0,hspace=0.3)
    fig.set_size_inches(6*ncols,5*nrows)
    plotname = "../plots/starrelations_"+suffix+"_"+plotname+".pdf"
    fig.savefig(plotname,bbox_inches='tight',dpi=150)

if __name__=="__main__":
    run(imfsims,"imf",["tracklength"],"sfe","tracklength")
    xvalues = ["compactness","firstmass","alltimemax","nphotons","nphotonstot","nphotonstff"]
    run(imfsims,"imf",xvalues,"sfe","all")

