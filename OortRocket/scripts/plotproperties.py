'''
Plot various properties over time in the simulations
Sam Geen, December 2017
'''

from startup import *

import inspect

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import matplotlib.patheffects as pe

import findproperties, starrelations, findHIIregions

from scipy.interpolate import interp1d

from matplotlib.lines import Line2D

import rdmfile

def _tsfeinsnap(snap):
    mgas = 1e4
    sink = sinks.FindSinks(snap)
    return np.sum(sink.mass) / mgas

tsfeinsnap = Hamu.Algorithm(_tsfeinsnap)

nprocs = 1

def tsfe(simname):
    global mgas
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, dum = starrelations.runforsim(simname,"firsttime")
    t, sfe = timefuncs.timefunc(sim,tsfeinsnap,processes=nprocs)
    if "MASS" in simname:
        sfe /= 10.0 # mgas is 10x larger here
    return t, sfe

momentuminsnap = Hamu.Algorithm(findproperties.totalmomentuminsnapS) # with sink momentum
def momentum(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t, p = timefuncs.timefunc(sim,momentuminsnap,processes=nprocs)
    #t -= tcreated
    return t, p

momentumatstarinsnap = Hamu.Algorithm(findproperties.radialmomentumatstarinsnap)
def momentumatstarpos(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t, p = timefuncs.timefunc(sim,momentumatstarinsnap,processes=nprocs)
    return t, p

massinsnap = Hamu.Algorithm(findproperties.massinsnap)
def mass(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    t, p = timefuncs.timefunc(sim,massinsnap,processes=nprocs)
    return t, p


neutralmassinsnap = Hamu.Algorithm(findproperties.neutralmassinsnap)
def neutralmass(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    t, m = timefuncs.timefunc(sim,neutralmassinsnap,processes=nprocs)
    return t, m

radiusinsnap = Hamu.Algorithm(findproperties.radiusinsnap3)
def radius(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,radiusinsnap,processes=nprocs)
    # HACK - fix boxlen error
    #boxlen = 0.121622418993404E+03
    #r *= boxlen
    #t -= tcreated
    return t, r

def windradiusratio(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,ri = timefuncs.timefunc(sim,radiusinsnap,processes=nprocs)
    t,rw = timefuncs.timefunc(sim,windradiusinsnap,processes=nprocs)
    # HACK - fix boxlen error
    #boxlen = 0.121622418993404E+03
    #r *= boxlen
    #t -= tcreated
    return t, rw / ri
    
windenergyinsnap = Hamu.Algorithm(findproperties.windenergyinsnap2)
def windenergy(simname,kinonly=False,thermonly=False):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t,e = timefuncs.timefunc(sim,windenergyinsnap,processes=1) # doesn't like tuples in the pool map
    #t -= tcreated
    therm = e[:,0]
    kin = e[:,1]
    tot = e[:,2]
    eout = np.array([tot, therm, kin]).T
    return t, eout

def windLemittedvscool(simname):
    fname = inspect.currentframe().f_code.co_name
    print "Running",fname,"for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t,Lwe = timefuncs.timefunc(sim,windLemittedHamu)
    t,Lwc = timefuncs.timefunc(sim,windLcoolinsnap,processes=nprocs)
    Lout = np.array([Lwe, Lwc]).T
    return t, Lout

windradiusinsnap = Hamu.Algorithm(findproperties.windradiusinsnap)
def windradius(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    if not "WIND" in simname:
        # Run some fast thing just to get zeros for each time
        t,dum = timefuncs.timefunc(sim,tsfeinsnap,processes=nprocs)
        return t,dum*0.0
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,windradiusinsnap,processes=nprocs)
    #t -= tcreated
    return t, r

windLcoolinsnap = Hamu.Algorithm(findproperties.windLcoolinsnap)
def windLcool(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,Lwc = timefuncs.timefunc(sim,windLcoolinsnap,processes=nprocs)
    #t -= tcreated
    return t, Lwc
    
radiusinsnapDEPRECATED = Hamu.Algorithm(findproperties.radiusinsnap)
def radiusDEPRECATED(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,radiusinsnapDEPRECATED)
    #t -= tcreated
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


def _windenergyemitted(snap):
    return starrelations.findcumulwinds(snap.hamusnap)[0]
windenergyemittedHamu = Hamu.Algorithm(_windenergyemitted)

def windenergyemitted(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, we = timefuncs.timefunc(sim,windenergyemittedHamu)
    return t, we

def _windLemitted(snap):
    return starrelations.findLwinds(snap.hamusnap)[0]
windLemittedHamu = Hamu.Algorithm(_windLemitted)

def windLemitted(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, we = timefuncs.timefunc(sim,windLemittedHamu)
    return t, we

def windenergyretained(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, wemit = timefuncs.timefunc(sim,windenergyemittedHamu)
    t, weingas = timefuncs.timefunc(sim,windenergyinsnap,processes=1) # doesn't like tuples in the pool map
    return t, weingas[:,2] / wemit

def _windmassemitted(snap):
    return starrelations.findcumulwinds(snap.hamusnap)[1]
windmassemittedHamu = Hamu.Algorithm(_windmassemitted)

def windmassemitted(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, wm = timefuncs.timefunc(sim,windmassemittedHamu)
    return t, wm

Pnontherminsnap = Hamu.Algorithm(findproperties.Pnontherminsnap)
def Pnontherm(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    t,r = timefuncs.timefunc(sim,Pnontherminsnap,verbose=True,processes=nprocs)
    return t, r

maxBfieldinsnap = Hamu.Algorithm(findproperties.maxBfieldinsnap)
def maxBfield(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    t,B = timefuncs.timefunc(sim,maxBfieldinsnap,verbose=True,processes=nprocs)
    return t, B

maxdensityinsnap = Hamu.Algorithm(findproperties.maxdensityinsnap)
def maxdensity(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    t,d = timefuncs.timefunc(sim,maxdensityinsnap,verbose=True,processes=nprocs)
    return t, d
#maxdensityinsnap.ForceReplaceCache(True) # Fucked up the last run

FindHIIRegions = findHIIregions.FindHIIRegions
def numHIIregions(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    t,n = timefuncs.timefunc(sim,FindHIIRegions,verbose=True,processes=nprocs)
    return t, n

photodensinsnap = Hamu.Algorithm(findproperties.photodensinsnap)
def photodens(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    t,pdens = timefuncs.timefunc(sim,photodensinsnap,verbose=True,processes=nprocs)
    return t, pdens

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

def run(simfunc,simnamesets,plotlabels,compare=False,secondfunc=None):
    plt.clf()
    numcols = len(simnamesets)
    wcloud = 1.71 # cloud density profile power law index
    funcname = simfunc.__name__
    fig, axes = plt.subplots(1,numcols,sharex=False,sharey=True)
    first = True
    rdm = rdmfile.RDMFile(__file__)
    for ax, simnames, plotlabel in zip(axes,simnamesets,plotlabels):
        linenames = []
        # Do func-related stuff for all plots
        ax.set_xscale("linear")
        ax.set_yscale("log")
        tlim = 0.1
        if not "MASS" in simnames[0]:
            tlim = 2.0
        #ax.set_xlabel("Time / Myr")
        ax.set_xlabel("Time after 1st star formed / Myr")
        textloc = 0.95
        xticks = None
        if funcname == "maxdensity":
            ax.set_ylim([0.1,1e11])
        if funcname == "nphotonsHII":
            #ax.set_xlim([1.5,12.0])
            ax.set_ylim([7e46,7e50])
            #xticks = [2,3,4,5,6,7,8,9,10]
        if funcname == "radius" or funcname == "radius2":
            #tlim = 1e-6
            if not compare:
                ax.set_ylim([0.05,100])
            #ax.set_xlabel("Time after 1st star formed / Myr")
            textloc = 0.1
        #if funcname == "windradius":
            #tlim = 1e-6
            #ax.set_xlabel("Time after 1st star formed / Myr")
        if funcname == "windenergy":
            #tlim = 1e-6
            ax.set_ylim([1e44,1e50])
            linenames = ["Total","Thermal","Kinetic"]
        if funcname == "windLemittedvscool":
            linenames = ["Emitted","Cooling"]
            ax.set_ylim([1e33,3e38])
        #if funcname == "momentum":
        #    tlim = 1e-6
        #    ax.set_xlabel("Time after 1st star formed / Myr")
        if funcname == "tsfe":
            ax.set_ylim([0.01,1.5])
            #xticks = [1,2,3,4,5,6,7,8,9,10]
            #ax.set_xlim([1,12.0])
        #if xticks is not None:
        #    xticklabels = [str(x) for x in xticks]
        #    ax.xaxis.set_major_locator(plt.NullLocator()) 
        #    ax.xaxis.set_minor_locator(plt.NullLocator()) 
        #    ax.set_xticks(xticks)
        #    ax.set_xticklabels(xticklabels)
        # Do main plot
        scale = ax.get_yscale()
        #if scale == "log":
        #    if funcname != "radius" and funcname != "momentum" and funcname != "radius2" and funcname != "windradius":
        #        ax.get_xaxis().set_minor_formatter(plt.ScalarFormatter()) 
        #        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter()) 
        if not compare:
            for simname in simnames:
                t, y = simfunc(simname)
                tcreated, sfe = starrelations.runforsim(simname,"firsttime")
                t -= tcreated
                #if scale == "log":
                #    mask = t > tlim
                #    t = t[mask]
                #    y = y[mask]
                legend2 = None
                if len(y.shape) == 1:
                    label = linestyles.Label(simname)
                    ax.plot(t,y,color=linestyles.Colour(simname),label=label,
                            linestyle=linestyles.Linestyle(simname),alpha=0.9,
                            path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
                    rdm.AddPoints(t,y,label=label)
                    if secondfunc:
                        t2, y2 = secondfunc(simname)
                        t2 -= tcreated
                        ax.plot(t2,y2,color=linestyles.Colour(simname),
                                linestyle=linestyles.Linestyle(simname),alpha=0.9,
                                path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()])
                        legelements = [Line2D([0],[0],color=linestyles.Colour(simname),
                                          linestyle="-",label=n,
                                          alpha=0.9,
                                          path_effects=[pe.Stroke(linewidth=lw, foreground='k'),
                                                        pe.Normal()]) for lw, n in zip([5,3], ["H II Region","Wind Bubble"])]
                        legend2 = ax.legend(handles=legelements, loc='upper right',framealpha=0.0,fontsize="x-small")
                else:
                    ntimes, nvals = y.shape
                    lines = ["-","--",":",":-"]
                    #import pdb; pdb.set_trace()
                    firstline = True
                    lines = lines[0:nvals]
                    legelements = [Line2D([0],[0],color=linestyles.Colour(simname),
                                          linestyle=l,label=n,
                                          alpha=0.9,
                                          path_effects=[pe.Stroke(linewidth=5, foreground='k'),
                                                        pe.Normal()]) for l, n in zip(lines, linenames)]
                    legend2 = ax.legend(handles=legelements, loc='lower center',framealpha=0.0)
                    for i in range(0,nvals):
                        if firstline:
                            label = linestyles.Label(simname)
                            firstline = False
                        else:
                            label = None
                        ax.plot(t,y[:,i],color=linestyles.Colour(simname),
                                label=label,
                                linestyle=lines[i],alpha=0.9,
                                path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
                        rdm.AddPoints(t,y[:,i],label=label)
        else:
            t1,y1 = simfunc(simnames[0])
            t2,y2 = simfunc(simnames[1])
            f = interp1d(t1,y1)
            tc = t2+0.0
            yc = y2 - f(tc)
            #if scale == "log":
            #    mask = tc > tlim
            #    tc = tc[mask]
            #    yc = yc[mask]
            tcreated, sfe = starrelations.runforsim(simnames[0],"firsttime")
            tc -= tcreated
            # +ve values (2nd = larger)
            names = [maxdensity,tsfe,mass,maxBfield,nphotonsHII,momentum,radius,
                     momentumatstarpos,windenergy,windradius,photodens,windLcool,windLemitted]
            effects = ["Maximum density","SFE","mass","max B field","Nphotons","momentum","radius",
                       "outflow momentum","wind bubble energy","wind bubble radius","average density of photoionised gas",
                       "Wind Bubble Cooling Luminosity","Wind Luminosity"]
            effectdict = {k:v for k,v in zip(names, effects)}
            effect = effectdict[simfunc]
            labelplus = "Winds increase "+effect
            labelminus = "Winds decrease "+effect
            ax.plot(tc,yc,color=linestyles.Colour(simnames[1]),label=labelplus,
                        path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
            # -ve values (2nd = smaller)
            ax.plot(tc,-yc,color=linestyles.Colour(simnames[1]),label=labelminus,
                    linestyle="--",
                    path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
            rdm.AddPoints(tc,yc,label=labelplus)
        # Overplot theoretical fits
        #if funcname == "momentum" and not compare:
            # Flat density profile
            #plotpowerlaw(ax,9.0/7.0,2e43,"k:",t0=0.5)
            # Power law density profile
            #plotpowerlaw(ax,(9.0-2.0*wcloud)/(7.0-2.0*wcloud),2e43,"k--",t0=0.5)
        #if (funcname == "radius" or funcname == "radius2") and not compare:
            # Flat density profile
            #plotpowerlaw(ax,4.0/7.0,10.0,"k:")
            # Power law density profile
            #plotpowerlaw(ax,4.0/(7.0-2.0*wcloud),10.0,"k--")
        # Set labels etc
        #ax.set_xlim([0.0,1])
        if funcname == "momentumatstarpos" or funcname == "radius" or funcname == "windradius":
            ax.set_xscale("log")
            #ax.set_xlim([3e-2,1])
        #if not "MASS" in simnames[0]:
        #    ax.set_xlim([3,7.3])
        #else:
        #    ax.set_xlim([1,5])
        ncollegend = 1
        if first:
            first = False
            legendloc = "upper left"
            # Dummy label to replace for final plots
            ax.set_ylabel(funcname)
            if funcname == "tsfe":
                ax.set_ylabel("SFE $\equiv M_{stars} / M_{ini}$")
            if funcname == "momentum":
                ax.set_ylabel("Momentum / g cm s$^{-1}$")
            if funcname == "momentumatstarpos":
                ax.set_ylabel("Outflow Momentum / g cm s$^{-1}$")
            if funcname == "radius":
                ax.set_ylabel("Mean HII Region radius / pc")
                legendloc = "upper left"
                ncollegend=3
            if funcname == "radius2":
                ax.set_ylabel("Sphericised HII Region radius / pc")
            if funcname == "nphotonsHII":
                ax.set_ylabel("Ionising Photon Emission Rate / s$^{-1}$")
            if funcname == "windradius":
                ax.set_ylabel("Sphericised radius of wind bubble / pc")
            if funcname == "windenergy":
                ax.set_ylabel("Energy in wind-shocked gas / erg")
                #legendloc = "lower right"
            if funcname == "windenergyemitted":
                ax.set_ylabel("Cumulative energy emitted by winds / erg")
            if funcname == "windLcool":
                ax.set_ylabel("Total cooling rate in wind bubble / erg/s")
            if funcname == "windLemitted":
                ax.set_ylabel("Wind luminosity / erg/s")
            if funcname == "windLemittedvscool":
                ax.set_ylabel("Wind luminosity (emitted vs cooling) / erg/s")
            if funcname == "windenergyretained":
                ax.set_ylabel("Fraction of wind energy input\nretained in the wind bubble",
                              multialignment='center')
                legendloc="lower right"
        if legendloc is not None:
            legend1 = ax.legend(fontsize="x-small",loc=legendloc,framealpha=0.0,ncol=ncollegend)
            legend1.get_frame().set_linewidth(0.0)
        if legend2 is not None:
            if funcname != "radius" or "DENSE" in simname:
                ax.add_artist(legend2)
                legend2.get_frame().set_linewidth(0.0)
        #ax.set_ylim([-0.01,0.25])
        ax.text(0.95, textloc,plotlabel, ha='right', va="top", transform=ax.transAxes)
    comparetxt = ""
    if compare:
        comparetxt = "_compare"
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(7*numcols,6)
    figname = plotfolder+funcname+"_both"+comparetxt+".pdf"
    fig.savefig(figname, dpi=80)
    rdm.Write(figname)
        
if __name__=="__main__":
    for func in [neutralmass]:
        run(func,(imfsims,icsims),(linestyles.starsc,linestyles.turbsc))
