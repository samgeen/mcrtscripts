'''
Plot various properties over time in the simulations
Sam Geen, December 2017
'''

from startup import *

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import matplotlib.patheffects as pe

import findproperties, starrelations, findHIIregions

from scipy.interpolate import interp1d

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

windradiusinsnap = Hamu.Algorithm(findproperties.windradiusinsnap)
def windradius(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
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

def run(simfunc,simnamesets,plotlabels,compare=False):
    plt.clf()
    numcols = len(simnamesets)
    wcloud = 1.71 # cloud density profile power law index
    funcname = simfunc.__name__
    fig, axes = plt.subplots(1,numcols,sharex=False,sharey=True)
    first = True
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
        if funcname == "nphotonsHII":
            #ax.set_xlim([1.5,12.0])
            ax.set_ylim([7e46,7e50])
            #xticks = [2,3,4,5,6,7,8,9,10]
        if funcname == "radius" or funcname == "radius2":
            #tlim = 1e-6
            if not compare:
                ax.set_ylim([0.1,70])
            #ax.set_xlabel("Time after 1st star formed / Myr")
            textloc = 0.1
        #if funcname == "windradius":
            #tlim = 1e-6
            #ax.set_xlabel("Time after 1st star formed / Myr")
        if funcname == "windenergy":
            #tlim = 1e-6
            #ax.set_ylim([0.99e45,1.1e49])
            linenames = ["Total","Thermal","Kinetic"]
        #if funcname == "momentum":
        #    tlim = 1e-6
        #    ax.set_xlabel("Time after 1st star formed / Myr")
        if funcname == "tsfe":
            ax.set_ylim([0.001,1.5])
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
                if len(y.shape) == 1:
                    ax.plot(t,y,color=linestyles.Colour(simname),label=linestyles.Label(simname),
                            linestyle=linestyles.Linestyle(simname),alpha=0.9,
                            path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
                else:
                    ntimes, nvals = y.shape
                    lines = ["-","--",":",":-"]
                    #import pdb; pdb.set_trace()
                    for i in range(0,nvals):
                        ax.plot(t,y[:,i],color=linestyles.Colour(simname),
                                label=linestyles.Label(simname)+" "+linenames[i],
                                linestyle=lines[i],alpha=0.9,
                                path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
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
            names = [tsfe,mass,maxBfield,nphotonsHII,momentum,radius,
                     momentumatstarpos,windenergy,windradius,photodens,windLcool]
            effects = ["SFE","mass","max B field","Nphotons","momentum","radius",
                       "outflow momentum","wind bubble energy","wind bubble radius","average density of photoionised gas",
                       "Wind Bubble Cooling Luminosity"]
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
        ax.set_xlim([-0.2,1])
        #if not "MASS" in simnames[0]:
        #    ax.set_xlim([3,7.3])
        #else:
        #    ax.set_xlim([1,5])
        if first:
            first = False
            legendloc = "upper left"
            # Dummy label to replace for final plots
            ax.set_ylabel(funcname)
            if funcname == "tsfe":
                ax.set_ylabel("SFE $\equiv M_{stars} / M_{ini}$")
            if funcname == "momentum":
                ax.set_ylabel("Momentum / g cm s$^{-1}$")
            if funcname == "radius":
                ax.set_ylabel("Mean HII Region radius / pc")
                legendloc = "center right"
            if funcname == "radius2":
                ax.set_ylabel("Sphericised HII Region radius / pc")
            if funcname == "nphotonsHII":
                ax.set_ylabel("Ionising Photon Emission Rate / s$^{-1}$")
            if funcname == "windradius":
                ax.set_ylabel("Sphericised radius of wind bubble / pc")
            if funcname == "windenergy":
                ax.set_ylabel("Total energy in wind-shocked gas / ergs")
        leg = ax.legend(fontsize="x-small",loc=legendloc)
        leg.get_frame().set_linewidth(0.0)
        #ax.set_ylim([-0.01,0.25])
        ax.text(0.95, textloc,plotlabel, ha='right', va="top", transform=ax.transAxes)
    comparetxt = ""
    if compare:
        comparetxt = "_compare"
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(7*numcols,6)
    fig.savefig(plotfolder+funcname+"_both"+comparetxt+".pdf", dpi=80)
        
if __name__=="__main__":
    #for func in [surfdens]:
    for func in [windLcool]:
        run(func,(["UVWIND_120","UVWIND_60","UVWIND_30"],["UVWIND_120_DENSE"]),
            ("Diffuse Cloud","Dense Cloud"))
    for compare in [True]:
        for func in [momentumatstarpos,tsfe,momentum,radius,nphotonsHII,photodens][::-1]:
            run(func,(["UV_120","UVWIND_120"],
                      ["UV_60","UVWIND_60"],
                      ["UV_30","UVWIND_30"],
                      ["UV_120_DENSE","UVWIND_120_DENSE"]),
                ("120 "+Msolar+" Star,\n Diffuse Cloud",
                 "60 "+Msolar+" Star,\n Diffuse Cloud",
                 "30 "+Msolar+" Star,\n Diffuse Cloud",
                 "120 "+Msolar+" Star,\n Dense Cloud"),compare=compare)
    for compare in [False]:
        for func in [momentumatstarpos,tsfe,momentum,radius,nphotonsHII,photodens]:
            run(func,(["NOFB","UV_120","UVWIND_120"],
                      ["NOFB","UV_60","UVWIND_60"],
                      ["NOFB","UV_30","UVWIND_30"],
                      ["NOFB_DENSE","UV_120_DENSE","UVWIND_120_DENSE"]),
                ("Diffuse Cloud","Diffuse Cloud","Diffuse Cloud","Dense Cloud"),compare=compare)
    for func in [windradiusratio,windenergyemitted,windmassemitted,windenergyretained,windenergy,windradius]:
        run(func,(["UVWIND_120","UVWIND_60","UVWIND_30"],["UVWIND_120_DENSE"]),
            ("Diffuse Cloud","Dense Cloud"))
