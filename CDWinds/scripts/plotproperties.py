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

sys.path.append("/home/stgeen0/MCRT/mcrtscripts/WindInUV/")

import solvecriteria
#import stars

Hamu.IGNOREERRORS = False
Hamu.ERRORCODE = -100.0

def _tsfeinsnap(snap):
    mgas = 1e4
    sink = sinks.FindSinks(snap)
    return np.sum(sink.mass) / mgas

tsfeinsnap = Hamu.Algorithm(_tsfeinsnap)

nprocs = 24

def tsfe(simname):
    global mgas
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, dum = starrelations.runforsim(simname,"firsttime")
    t, sfe = timefuncs.timefunc(sim,tsfeinsnap,processes=nprocs)
    if "MASS" in simname:
        sfe /= 10.0 # mgas is 10x larger here
    return t, sfe

momentuminsnap = Hamu.Algorithm(findproperties.totalmomentuminsnapS) # with sink momentum
def momentum(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t, p = timefuncs.timefunc(sim,momentuminsnap,processes=nprocs)
    #t -= tcreated
    return t, p

momentumatstarinsnap = Hamu.Algorithm(findproperties.radialmomentumatstarinsnap)
def momentumatstarpos(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t, p = timefuncs.timefunc(sim,momentumatstarinsnap,processes=nprocs)
    return t, p

massinsnap = Hamu.Algorithm(findproperties.massinsnap)
def mass(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    t, p = timefuncs.timefunc(sim,massinsnap,processes=nprocs)
    return t, p
    
radiusinsnap = Hamu.Algorithm(findproperties.radiusinsnap3)
def radius(sim):
    if type(sim) == type("duck"):
        sim = hamusims[sim]
    simname = sim.Name()
    print("Running for simulation", simname)
    #tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,radiusinsnap,processes=nprocs)
    # HACK - fix boxlen error
    #boxlen = 0.121622418993404E+03
    #r *= boxlen
    #t -= tcreated
    return t, r

maxradiusatstarposinsnap = Hamu.Algorithm(findproperties.maxradiusatstarpos)
def maxradiusatstarpos(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,maxradiusatstarposinsnap,processes=nprocs)
    return t, r

maxwindradiusatstarposinsnap = Hamu.Algorithm(findproperties.maxwindradiusatstarpos)
def maxwindradiusatstarpos(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,maxwindradiusatstarposinsnap,processes=nprocs)
    return t, r


energyinsnap = Hamu.Algorithm(findproperties.energyinsnap)
def energy(simname,kinonly=False,thermonly=False):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t,e = timefuncs.timefunc(sim,energyinsnap,processes=1) # doesn't like tuples in the pool map
    #t -= tcreated
    therm = e[:,0]
    kin = e[:,1]
    tot = e[:,2]
    eout = np.array([tot, therm, kin]).T
    return t, eout

windpressureinsnap = Hamu.Algorithm(findproperties.windpressureinsnap)
def windpressure(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,windpressureinsnap,processes=nprocs)
    return t, r

def windradiusratio(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,ri = timefuncs.timefunc(sim,radiusinsnap,processes=nprocs)
    # Don't waste time on this if there are no winds in the run
    if "WIND" in simname:
        t,rw = timefuncs.timefunc(sim,windradiusinsnap,processes=nprocs)
    else:
        rw = ri * 0.0
    # HACK - fix boxlen error
    #boxlen = 0.121622418993404E+03
    #r *= boxlen
    #t -= tcreated
    return t, rw / ri

def windradiusratio_analytic(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    metal = 0.014
    starmass = float(str.join("",[i for i in simname if i.isdigit()]))
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    ts,ris = timefuncs.timefunc(sim,radiusinsnap,processes=nprocs)
    star = stars.Star(starmass,metal)
    rwvsri = []
    for t, ri in zip(ts, ris):
        ci = solvecriteria.SoundSpeedFromT(1e4) # from sampling simulations, is very close to this
        print("ci", ci)
        cloud = solvecriteria.CloudFromr0n0(ri*pcincm,432.0,metal,ci=ci) # density doesn't really matter for this
        tins = (t-tcreated)*Myrins
        rivsrw = solvecriteria.RivsRw(tins,star,cloud)
        if rivsrw > 0.0:
            rwvsri.append(1.0/rivsrw)
        else:
            rwvsri.append(0.0)
    return ts, np.array(rwvsri)
    
windenergyinsnap = Hamu.Algorithm(findproperties.windenergyinsnap3)
def windenergy(sim,kinonly=False,thermonly=False):
    if type(sim) == type("duck"):
        sim = hamusims[sim]
    simname = sim.Name()
    print("Running for simulation", simname)
    #tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t,e = timefuncs.timefunc(sim,windenergyinsnap,processes=1) # doesn't like tuples in the pool map
    print(e.shape)
    #t -= tcreated
    therm = e[:,0]
    kin = e[:,1]
    tot = e[:,2]
    eout = np.array([tot, therm, kin]).T
    return t, eout

def windLemittedvscool(simname):
    fname = inspect.currentframe().f_code.co_name
    print("Running",fname,"for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t,Lwe = timefuncs.timefunc(sim,windLemittedHamu)
    t,Lwc = timefuncs.timefunc(sim,windLcoolinsnap,processes=nprocs)
    t,Lwx = timefuncs.timefunc(sim,xrayLcoolinsnap,processes=nprocs)
    Lout = np.array([Lwe, Lwc, Lwx]).T
    return t, Lout

windradiusinsnap = Hamu.Algorithm(findproperties.windradiusinsnap)
def windradius(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    if not "WIND" in simname:
        # Run some fast thing just to get zeros for each time
        t,dum = timefuncs.timefunc(sim,tsfeinsnap,processes=nprocs)
        return t,dum*0.0
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,windradiusinsnap,processes=nprocs)
    #t -= tcreated
    return t, r

freestreamradiusinsnap = Hamu.Algorithm(findproperties.freestreamradiusinsnap)
def freestreamradius(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    if not "WIND" in simname:
        # Run some fast thing just to get zeros for each time
        t,dum = timefuncs.timefunc(sim,tsfeinsnap,processes=nprocs)
        return t,dum*0.0
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,freestreamradiusinsnap,processes=nprocs)
    #t -= tcreated
    return t, r

windLcoolinsnap = Hamu.Algorithm(findproperties.windLcoolinsnap)
xrayLcoolinsnap = Hamu.Algorithm(findproperties.xrayLcoolinsnap)
def windLcool(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,Lwc = timefuncs.timefunc(sim,windLcoolinsnap,processes=nprocs)
    #t -= tcreated
    return t, Lwc
    
radiusinsnapDEPRECATED = Hamu.Algorithm(findproperties.radiusinsnap)
def radiusDEPRECATED(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")  
    t,r = timefuncs.timefunc(sim,radiusinsnapDEPRECATED)
    #t -= tcreated
    return t, r

def _NphotonsHII(snap):
    return np.sum(starrelations.findNphotons(snap.hamusnap)[3:5])
NphotonsHII = Hamu.Algorithm(_NphotonsHII)

def nphotonsHII(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, nphot = timefuncs.timefunc(sim,NphotonsHII)
    #t -= tcreated
    return t, nphot


def _windenergyemitted(snap):
    return starrelations.findcumulwinds(snap.hamusnap)[0]
windenergyemittedHamu = Hamu.Algorithm(_windenergyemitted)

def windenergyemitted(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, we = timefuncs.timefunc(sim,windenergyemittedHamu)
    return t, we

def _windLemitted(snap):
    return starrelations.findLwinds(snap.hamusnap)[0]
windLemittedHamu = Hamu.Algorithm(_windLemitted)

def windLemitted(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, we = timefuncs.timefunc(sim,windLemittedHamu)
    return t, we


def _windmomemitted(snap):
    return starrelations.findcumulmomwinds(snap.hamusnap)
windmomemittedHamu = Hamu.Algorithm(_windmomemitted)

def windmomemitted(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, wm = timefuncs.timefunc(sim,windmomemittedHamu,processes=nprocs)
    #import pdb; pdb.set_trace()
    return t, wm


def windenergyretained(simname,sim=None):
    print("Running for simulation", simname)
    if sim is None:
        sim = hamusims[simname]
    #tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    #tcreated = timefuncs.FindTcreatedFirstStar(sim)
    t, wemit = timefuncs.timefunc(sim,windenergyemittedHamu)
    t, weingas = timefuncs.timefunc(sim,windenergyinsnap,processes=1) # doesn't like tuples in the pool map
    print(weingas, weingas.shape, wemit, t)
    return t, weingas[:,2] / wemit

def _windmassemitted(snap):
    return starrelations.findcumulwinds(snap.hamusnap)[1]
windmassemittedHamu = Hamu.Algorithm(_windmassemitted)

def windmassemitted(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t, wm = timefuncs.timefunc(sim,windmassemittedHamu)
    return t, wm

Pnontherminsnap = Hamu.Algorithm(findproperties.Pnontherminsnap)
def Pnontherm(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    t,r = timefuncs.timefunc(sim,Pnontherminsnap,verbose=True,processes=nprocs)
    return t, r

maxBfieldinsnap = Hamu.Algorithm(findproperties.maxBfieldinsnap)
def maxBfield(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    t,B = timefuncs.timefunc(sim,maxBfieldinsnap,verbose=True,processes=nprocs)
    return t, B

maxdensityinsnap = Hamu.Algorithm(findproperties.maxdensityinsnap)
def maxdensity(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    t,d = timefuncs.timefunc(sim,maxdensityinsnap,verbose=True,processes=nprocs)
    return t, d
#maxdensityinsnap.ForceReplaceCache(True) # Fucked up the last run

FindHIIRegions = findHIIregions.FindHIIRegions
def numHIIregions(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    t,n = timefuncs.timefunc(sim,FindHIIRegions,verbose=True,processes=nprocs)
    return t, n

photodensinsnap = Hamu.Algorithm(findproperties.photodensinsnap)
def photodens(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    t,pdens = timefuncs.timefunc(sim,photodensinsnap,verbose=True,processes=nprocs)
    return t, pdens

Bfieldenergyinsnap = Hamu.Algorithm(findproperties.Bfieldenergyinsnap2)
def Bfieldenergy(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    t, BE = timefuncs.timefunc(sim,Bfieldenergyinsnap,processes=nprocs)
    # Correct BE since we use microGauss as internal units, not Gauss
    return t, BE

def energyplusB(simname):
    print("Running for simulation", simname)
    sim = hamusims[simname]
    tcreated, sfe = starrelations.runforsim(simname,"firsttime")
    t,e = timefuncs.timefunc(sim,energyinsnap,processes=1) # doesn't like tuples in the pool map
    t,BE = timefuncs.timefunc(sim,Bfieldenergyinsnap,processes=1) # doesn't like tuples in the pool map
    #t -= tcreated
    therm = e[:,0]
    kin = e[:,1]
    eout = np.array([BE, therm, kin]).T
    return t, eout

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

def run(simfunc,simnamesets,plotlabels,compare=False,secondfuncs=None,gradient=False,suffix=""):
    plt.clf()
    linestyles.reset()
    numcols = len(simnamesets)
    wcloud = 1.71 # cloud density profile power law index
    funcname = simfunc.__name__
    fig, axes = plt.subplots(1,numcols,sharex=False,sharey=True)
    if len(simnamesets) == 1:
        axes = [axes,]
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
        ax.set_xlabel("Time after star formed / Myr")
        textloc = 0.95
        xticks = None
        if funcname == "maxdensity":
            ax.set_ylim([0.1,1e11])
        if funcname == "nphotonsHII":
            #ax.set_xlim([1.5,12.0])
            ax.set_ylim([7e46,7e50])
            #xticks = [2,3,4,5,6,7,8,9,10]
        if funcname == "momentumatstarpos" and not gradient:
            textloc = 0.1    
        if "radius" in funcname:
            #tlim = 1e-6
            #if not compare:
            #    ax.set_ylim([0.08,110])
            #ax.set_xlabel("Time after 1st star formed / Myr")
            textloc = 0.1
            ax.set_xscale("linear")
            ax.set_yscale("linear")
        # HARD-CODE TO ALLOW COMPARISON
        if "maxradius" in funcname:
            for setname in simnamesets:
                if "physics" in setname:
                    ax.set_ylim([0,13])
                if "seeds" in setname:
                    ax.set_ylim([0,25])
        #if funcname == "windradius":
            #tlim = 1e-6
            #ax.set_xlabel("Time after 1st star formed / Myr")
        if funcname == "windenergy":
            #tlim = 1e-6
            ax.set_ylim([1e44,1e50])
            linenames = ["Total","Thermal","Kinetic"]
        if funcname == "energy":
            linenames = ["Total","Thermal","Kinetic"]
            if funcname == "energyplusB":
                linenames = ["Magnetic","Thermal","Kinetic"]
        if funcname == "windLemittedvscool":
            linenames = ["Emitted","Cooling","Cooling $> 10^{6}$ K"]
            ax.set_ylim([1e31,1e39])
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
                    if secondfuncs:
                        secondlines = ["--",":","-."]
                        # HACK, just use the simname linestyles
                        secondlines = [linestyles.Linestyle(simname)]*3
                        for i, secondfunc in enumerate(secondfuncs):
                            secondline = secondlines[i]
                            t2, y2 = secondfunc(simname)
                            t2 -= tcreated
                            ax.plot(t2,y2,color=linestyles.Colour(simname),
                                    linestyle=secondline,alpha=0.9,
                                    path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()])
                        if funcname == "radius":
                            rlabels = ["H II Region","Wind Bubble"]#,"$V >$ 1000 km/s"]
                            legelements = [Line2D([0],[0],color=linestyles.Colour(simname),
                                                  linestyle=ls,label=n,
                                                  alpha=0.9,
                                                  path_effects=[pe.Stroke(linewidth=lw, foreground='k'),
                                                                pe.Normal()]) for ls,lw, n in zip(["-"]+secondlines,
                                                                                                  [5,3,3],
                                                                                                  rlabels)]
                            legend2 = ax.legend(handles=legelements, loc='upper right',framealpha=0.0,fontsize="x-small")
                        if funcname == "windradiusratio":
                            #import pdb; pdb.set_trace()
                            legelements = [Line2D([0],[0],color=linestyles.Colour(simname),
                                                  linestyle=ls,label=n,
                                                  alpha=0.9,
                                                  path_effects=[pe.Stroke(linewidth=lw, foreground='k'),
                                                                pe.Normal()]) for ls,lw, n in zip(["-",secondlines[0]],
                                                                                                  [5,3],
                                                                                                  ["Simulation",
                                                                                                   "Analytic Model"])]
                            legend2 = ax.legend(handles=legelements, loc='lower right',framealpha=0.0,fontsize="x-small")
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
                    # Lol hacks
                    if funcname != "windLemittedvscool" or "DENSE" in simname:
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
            legend2 = None
            t1,y1 = simfunc(simnames[0])
            t2,y2 = simfunc(simnames[1])
            f = interp1d(t1,y1)
            tc = t2+0.0
            yc = y2[tc > t1.min()]
            tc = tc[tc > t1.min()]
            yc = yc - f(tc)
            #if scale == "log":
            #    mask = tc > tlim
            #    tc = tc[mask]
            #    yc = yc[mask]
            tcreated, sfe = starrelations.runforsim(simnames[0],"firsttime")
            tc -= tcreated
            # +ve values (2nd = larger)
            names = [maxdensity,tsfe,mass,maxBfield,nphotonsHII,momentum,radius,
                     momentumatstarpos,windenergy,windradius,freestreamradius,photodens,windLcool,windLemitted,
                     windpressure]
            effects = ["Maximum density","SFE","mass","max B field","Nphotons","momentum","radius",
                       "outflow momentum","wind bubble energy","wind bubble radius","free streaming radius",
                       "average density of photoionised gas",
                       "Wind Bubble Cooling Luminosity","Wind Luminosity","Pressure in Wind Bubble"]
            effectdict = {k:v for k,v in zip(names, effects)}
            effect = effectdict[simfunc]
            labelplus = "Winds increase "+effect
            labelminus = "Winds decrease "+effect
            if gradient:
                newt = np.linspace(tc[0],tc[-1],1000)
                newy = interp1d(tc,yc)(newt)
                tc = newt
                yc = np.concatenate(([0],np.diff(newy)/np.diff(newt))).flatten()/Myrins
            ax.plot(tc,yc,color=linestyles.Colour(simnames[1]),label=labelplus,
                        path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
            # -ve values (2nd = smaller)
            ax.plot(tc,-yc,color=linestyles.Colour(simnames[1]),label=labelminus,
                    linestyle="--",
                    path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
            rdm.AddPoints(tc,yc,label=labelplus)
            if secondfuncs:
                secondlines = [":","-."]
                # HACK, just use the simname linestyles
                secondlines = [linestyles.Linestyle(simnames[1])]*2
                for i, secondfunc in enumerate(secondfuncs):
                    secondline = secondlines[i]
                    t2, y2 = secondfunc(simnames[1])
                    t2 -= tcreated
                    if gradient:
                        newt = np.linspace(t2[0],t2[-1],1000)
                        newy = interp1d(t2,y2)(newt)
                        t2 = newt
                        y2 = np.concatenate(([0],np.diff(newy)/np.diff(newt))).flatten()/Myrins
                    ax.plot(t2,y2,color=linestyles.Colour(simnames[1]),
                            linestyle=secondline,alpha=0.9,
                            path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()])
                    if gradient:
                        ax.set_ylim([4e40/Myrins,1.2e44/Myrins])
                    else:
                        ax.set_ylim([1e40,1e42])
                    if funcname == "momentumatstarpos":
                            rlabels = ["$\dot{p}_{\mathrm{w}}$"]#,"$V >$ 1000 km/s"]
                            legelements = [Line2D([0],[0],color=linestyles.Colour(simnames[1]),
                                                  linestyle=ls,label=n,
                                                  alpha=0.9) for ls,lw, n in zip([":"]+secondlines,
                                                                                 [5,3,3],
                                                                                 rlabels)]
                            legloc2 = "lower right"
                            if "30" in simnames[1]:
                                legloc2 = "center left"
                            legend2 = ax.legend(handles=legelements, loc=legloc2,framealpha=0.0,fontsize="small")
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
        ax.set_xlim([0.0,0.4])
        #if funcname == "momentumatstarpos" or funcname == "radius" or funcname == "windradius" or funcname == "freestreamradius":
            #ax.set_xscale("log")
            #ax.set_xlim([3e-2,0.4])
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
                if gradient:
                    ax.set_ylabel("$\mathrm{d}\Delta p_{(\mathrm{UVWIND-UV})} / \mathrm{dt}$ / g cm s$^{-2}$")
                else:
                    ax.set_ylabel("Outflow Momentum / g cm s$^{-1}$")
                    ncollegend = 2
            if funcname == "radius":
                ax.set_ylabel("Mean HII Region radius / pc")
                legendloc = "upper left"
                #ncollegend=3
            if funcname == "radius2":
                ax.set_ylabel("Sphericised HII Region radius / pc")
            if funcname == "nphotonsHII":
                ax.set_ylabel("Ionising Photon Emission Rate / s$^{-1}$")
            if funcname == "windradius":
                ax.set_ylabel("Sphericised radius of wind bubble / pc")
            if funcname == "windradiusratio":
                ax.set_ylabel("$r_{w,s} / r_{i,s}$")
                legendloc = "upper center"
            if funcname == "windenergy":
                ax.set_ylabel("Energy in wind-shocked gas / erg")
                #legendloc = "lower right"
            if funcname == "windenergyemitted":
                ax.set_ylabel("Cumulative energy emitted by winds / erg")
            if funcname == "windpressure":
                ax.set_ylabel("Pressure in wind bubble / erg / cm$^3$")
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
        comparetxt += "_compare"
    if gradient:
        comparetxt += "_gradient"
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(8*numcols,6)
    figname = plotfolder+funcname+"_both"+comparetxt+suffix+".pdf"
    fig.savefig(figname, dpi=80)
    rdm.Write(figname)
    fig.clear()
        
def runall():
    '''
    for compare in [True]:
        for func in [momentumatstarpos,tsfe,momentum,radius,nphotonsHII,photodens]:
            run(func,(["UV_30","UVWIND_30"],
                      ["UV_60","UVWIND_60"],
                      ["UV_120","UVWIND_120"],
                      ["UV_120_DENSE","UVWIND_120_DENSE"]),
                ("30 "+Msolar+" Star,\n Diffuse Cloud",
                 "60 "+Msolar+" Star,\n Diffuse Cloud",
                 "120 "+Msolar+" Star,\n Diffuse Cloud",
                 "120 "+Msolar+" Star,\n Dense Cloud"),compare=compare)
    '''
    # [momentumatstarpos,tsfe,momentum,radius,nphotonsHII,photodens]
    #alluvnames = ["UV_"+str(num) for num in [30,60,120]]
    #allwindnames = ["UVWIND_"+str(num) for num in [30,60,120]]
    #allwindpressnames = ["UVWINDPRESS_"+str(num) for num in [30,60,120]]
    #allfbnames = alluvnames+allwindnames+allwindpressnames
    #allnames = ["NOFB"]+allfbnames
    #windpressnames = ["UVWIND_120_DENSE"]

    #for func in [maxBfield]:
    #    run(func,(["NOFB"],["NOFB_DENSE"]),
    #        ("Diffuse Cloud","Dense Cloud"),compare=False)

    for setname, simset in simsets.items():
        linestyles.CURRSIMSET = setname
        #setname = "physics"
        #simset = physicsset
        #    if True:

        #for func in [maxradiusatstarpos,maxwindradiusatstarpos,
        #             windLemittedvscool,windenergyemitted,windmassemitted,
        #             windenergyretained,windenergy,windradius,freestreamradius]:
        #    run(func,[simset,],
        #        [setname,],compare=False,suffix=setname)

        #for func in [windenergyretained]:
            #,maxradiusatstarpos,maxwindradiusatstarpos,
            #         windenergyemitted,windmassemitted,
            #         windenergyretained,windenergy,windradius,freestreamradius]:

        for func in [windradiusratio,radius,windpressure,momentumatstarpos,maxradiusatstarpos,maxwindradiusatstarpos,
                     windLemittedvscool,windenergyemitted,windmassemitted,
                     windenergyretained,windenergy,windradius,freestreamradius]:
            run(func,[simset,],
                [setname,],compare=False,suffix=setname)

        for func in [energyplusB]:
            run(func,[simset,], # ,"UVWINDPRESS_120_DENSE"]),
                [setname,],compare=False,suffix=setname)

        for func in [Bfieldenergy]:
            run(func,[simset,], # ,"UVWINDPRESS_120_DENSE"]),
                [setname,],compare=False,suffix=setname)

    '''
        
        for func in [momentumatstarpos]:
            run(func,(["UV_30","UVWIND_30"],
                    ["UV_60","UVWIND_60"],
                    ["UV_120","UVWIND_120"],
                    ["UV_120_DENSE","UVWIND_120_DENSE"]),
                ("30 "+Msolar+" Star,\n Diffuse Cloud",
                "60 "+Msolar+" Star,\n Diffuse Cloud",
                "120 "+Msolar+" Star,\n Diffuse Cloud",
                "120 "+Msolar+" Star,\n Dense Cloud"),compare=True,secondfuncs=(windmomemitted,))
            run(func,(["UVWINDPRESS_120","UVWIND_120","UV_120",
                    "UVWINDPRESS_60","UVWIND_60","UV_60",
                    "UVWINDPRESS_30","UVWIND_30","UV_30",
                    "NOFB"],
                    ["NOFB_DENSE","UV_120_DENSE","UVWIND_120_DENSE","UVWINDPRESS_120_DENSE"][::-1]),
                ("Diffuse Cloud","Dense Cloud"),compare=False)
        for func in [momentumatstarpos]:
            run(func,(["UV_30","UVWIND_30"],
                    ["UV_60","UVWIND_60"],
                    ["UV_120","UVWIND_120"],
                    ["UV_120_DENSE","UVWIND_120_DENSE"]),
                ("30 "+Msolar+" Star,\n Diffuse Cloud",
                "60 "+Msolar+" Star,\n Diffuse Cloud",
                "120 "+Msolar+" Star,\n Diffuse Cloud",
                "120 "+Msolar+" Star,\n Dense Cloud"),compare=True,secondfuncs=(windmomemitted,))
            run(func,(["UV_30","UVWIND_30"],
                    ["UV_60","UVWIND_60"],
                    ["UV_120","UVWIND_120"],
                    ["UV_120_DENSE","UVWIND_120_DENSE"]),
                ("30 "+Msolar+" Star,\n Diffuse Cloud",
                "60 "+Msolar+" Star,\n Diffuse Cloud",
                "120 "+Msolar+" Star,\n Diffuse Cloud",
                "120 "+Msolar+" Star,\n Dense Cloud"),compare=True,secondfuncs=(windmomemitted,),gradient=True)
            # Same but split up
            run(func,(["UV_30","UVWIND_30"],
                    ["UV_60","UVWIND_60"]),
                ("30 "+Msolar+" Star,\n Diffuse Cloud",
                "60 "+Msolar+" Star,\n Diffuse Cloud"),
                compare=True,secondfuncs=(windmomemitted,),gradient=True,suffix="1of2")
            run(func,(["UV_120","UVWIND_120"],
                    ["UV_120_DENSE","UVWIND_120_DENSE"]),
                ("120 "+Msolar+" Star,\n Diffuse Cloud",
                "120 "+Msolar+" Star,\n Dense Cloud"),
                compare=True,secondfuncs=(windmomemitted,),gradient=True,suffix="2of2")
                    
        for func in [windradiusratio]:
            run(func,(allwindpressnames,["UVWINDPRESS_120_DENSE"]),
                ("Diffuse Cloud","Dense Cloud"),compare=False,secondfuncs=(windradiusratio_analytic,))
            

        for func in [radius]:
            sfuncs = (windradius,) #freestreamradius)
            run(func,(allfbnames,
                    ["UV_120_DENSE","UVWIND_120_DENSE","UVWINDPRESS_120_DENSE"]),
                ("Diffuse Cloud","Dense Cloud"),compare=False,secondfuncs=sfuncs)

        for func in [windLemittedvscool,windenergyemitted,windmassemitted,windenergyretained,windenergy,windradius,freestreamradius]:
            run(func,(allwindpressnames,["UVWINDPRESS_120_DENSE"]),
                ("Diffuse Cloud","Dense Cloud"))

        '''

if __name__=="__main__":
    runall()        
