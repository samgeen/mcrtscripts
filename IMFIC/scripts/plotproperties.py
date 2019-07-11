'''
Find the total star formation efficiency
Sam Geen, December 2017
'''

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

from startup import *

import matplotlib.patheffects as pe

import findproperties, starrelations

def _tsfeinsnap(snap):
    sink = sinks.FindSinks(snap)
    mgas = 1e4
    return np.sum(sink.mass) / mgas

tsfeinsnap = Hamu.Algorithm(_tsfeinsnap)

def tsfe(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t, sfe = timefuncs.timefunc(sim,tsfeinsnap)
    #t -= tcreated
    return t, sfe

momentuminsnap = Hamu.Algorithm(findproperties.totalmomentuminsnap)
def momentum(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t, p = timefuncs.timefunc(sim,momentuminsnap)
    t -= tcreated
    return t, p

radiusinsnap2 = Hamu.Algorithm(findproperties.radiusinsnap3)
def radius2(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,radiusinsnap2)
    # HACK - fix boxlen error
    #boxlen = 0.121622418993404E+03
    #r *= boxlen
    t -= tcreated
    return t, r

radiusinsnap = Hamu.Algorithm(findproperties.radiusinsnap)
def radius(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,radiusinsnap)
    t -= tcreated
    return t, r

def _NphotonsHII(snap):
    return np.sum(starrelations.findNphotons(snap.hamusnap)[3:5])
NphotonsHII = Hamu.Algorithm(_NphotonsHII)

def nphotonsHII(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t, nphot = timefuncs.timefunc(sim,NphotonsHII)
    #t -= tcreated
    return t, nphot
    

def plotpowerlaw(ax,w,y0,linestyle,t0=0.3):
    '''
    Plot a power law fit fit ( p = c * t^w )
    Manually pick a y intercept at y0 for t=3.0 Myr
    '''
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    t = np.arange(xlims[0],xlims[1],0.02)
    #y0 = 2e43
    #w = 9.0/7.0
    c = y0 / t0**w
    x = t
    y = c * t**w
    mask = (y > ylims[0])*(y < ylims[1])
    x = x[mask]
    y = y[mask]
    ax.plot(x,y,linestyle,zorder=0,alpha=0.7)

def run(simfunc,simnamesets,plotlabels):
    plt.clf()
    wcloud = 1.71 # cloud density profile power law index
    funcname = simfunc.__name__
    fig, axes = plt.subplots(1,2,sharex=True,sharey=True)
    first = True
    for ax, simnames, plotlabel in zip(axes,simnamesets,plotlabels):
        # Do func-related stuff for all plots
        ax.set_xscale("log")
        ax.set_yscale("log")
        tlim = 2.0
        ax.set_xlabel("Time / Myr")
        textloc = 0.95
        xticks = None
        if funcname == "nphotonsHII":
            ax.set_xlim([1.5,12.0])
            ax.set_ylim([7e45,7e49])
            xticks = [2,3,4,5,6,7,8,9,10]
        if funcname == "radius" or funcname == "radius2":
            tlim = 1e-6
            ax.set_ylim([0.1,70])
            ax.set_xlabel("Time after 1st star formed / Myr")
            textloc = 0.1
        if funcname == "momentum":
            tlim = 1e-6
            ax.set_xlabel("Time after 1st star formed / Myr")
        if funcname == "tsfe":
            xticks = [2,3,4,5,6,7,8,9,10]
            ax.set_xlim([1.5,12.0])
        if xticks is not None:
            xticklabels = [str(x) for x in xticks]
            ax.xaxis.set_major_locator(plt.NullLocator()) 
            ax.xaxis.set_minor_locator(plt.NullLocator()) 
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
        # Do main plot
        for simname in simnames:
            t, y = simfunc(simname)
            scale = ax.get_yscale()
            if scale == "log":
                if funcname != "radius" and funcname != "momentum" and funcname != "radius2":
                    ax.get_xaxis().set_minor_formatter(plt.ScalarFormatter()) 
                    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter()) 
                mask = t > tlim
                t = t[mask]
                y = y[mask]
            ax.plot(t,y,color=linestyles.Colour(simname),label=linestyles.Label(simname),
                    linestyle=linestyles.Linestyle(simname),
                    path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
        # Overplot theoretical fits
        if funcname == "momentum":
            # Flat density profile
            plotpowerlaw(ax,9.0/7.0,2e43,"k:",t0=0.5)
            # Power law density profile
            plotpowerlaw(ax,(9.0-2.0*wcloud)/(7.0-2.0*wcloud),2e43,"k--",t0=0.5)
        if funcname == "radius" or funcname == "radius2":
            # Flat density profile
            plotpowerlaw(ax,4.0/7.0,10.0,"k:")
            # Power law density profile
            plotpowerlaw(ax,4.0/(7.0-2.0*wcloud),10.0,"k--")
        # Set labels etc
        if first:
            if funcname == "tsfe":
                ax.set_ylabel("SFE $\equiv M_{stars} / M_{ini}$")
            if funcname == "momentum":
                ax.set_xlim([6e-3,12])
                ax.set_ylabel("Momentum / g cm s$^{-1}$")
            if funcname == "radius":
                ax.set_ylabel("Mean HII Region radius / pc")
            if funcname == "radius2":
                ax.set_ylabel("Sphericised HII Region radius / pc")
            if funcname == "nphotonsHII":
                ax.set_ylabel("Ionising Photon Emission Rate / s$^{-1}$")
        leg = ax.legend(fontsize="x-small",loc="upper left")
        leg.get_frame().set_linewidth(0.0)
        #ax.set_ylim([-0.01,0.25])
        ax.text(0.95, textloc,plotlabel, ha='right', va="top", transform=ax.transAxes)
        first = False
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(14,6)
    fig.savefig(plotfolder+funcname+"_both.pdf", dpi=80)
        
if __name__=="__main__":
    for func in [momentum,nphotonsHII,tsfe,radius,radius2]:
        run(func,(imfsims,icsims),(linestyles.starsc,linestyles.turbsc))
