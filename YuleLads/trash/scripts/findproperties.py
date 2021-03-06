'''
Plot properties in simulations
Sam Geen, November 2015
'''

import os, glob, time
from adjustText import adjust_text

from startup import *

import Hamu
import pymses
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from pymses.filters import CellsToPoints
from pymses.utils import constants as C

import linestyles
import sinks, stellars
import skimfunc

import copy

gamma = 1.4 # Ramses default
mH = 1.66e-24
kB = 1.38e-16

def stellarmassinsnap(snap):
    print "Finding stellar mass in snapshot"
    sink = sinks.FindSinks(snap.hamusnap)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    return np.sum(sink.mass)

def maxstellarmassinsnap(snap):
    print "Finding max stellar mass in snapshot"
    stellar = stellars.FindStellar(snap.hamusnap)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    if len(stellar.mass) > 0:
        mmax = np.max(stellar.mass)
    else:
        mmax = 0.0
    return mmax

def massinsnap(snap,nHthresh=0.0):
    print "Finding mass in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cells = CellsToPoints(amr).flatten()
    vols = (cells.get_sizes())**3.0
    dens = cells["rho"]
    mass = dens*vols*snap.info["unit_mass"].express(C.Msun)
    if nHthresh > 0.0:
        mass = mass[dens > nHthresh]
    return np.sum(mass)

def meanionpressinsnap(snap):
    print "Finding mean pressure of ionised gas in snap", snap.iout
    amr = snap.amr_source(["P","xHII"])
    cells = CellsToPoints(amr).flatten()
    vols = (cells.get_sizes())**3.0
    P = cells["P"]
    E = P*vols
    mask = cells["xHII"] > 0.5
    meanpress = np.sum(E[mask]) / np.sum(vols[mask])
    return meanpress

def meaniondensinsnap(snap):
    print "Finding mean density of ionised gas in snap", snap.iout
    amr = snap.amr_source(["rho","xHII"])
    cells = CellsToPoints(amr).flatten()
    vols = (cells.get_sizes())**3.0
    dens = cells["rho"]
    mass = dens*vols*snap.info["unit_mass"].express(C.Msun)
    mask = cells["xHII"] > 0.5
    meandens = np.sum(mass[mask]) / np.sum(vols[mask])
    return meandens

def etherminsnap(snap):
    print "Finding max T in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cells = CellsToPoints(amr).flatten()
    vols = (cells.get_sizes())**3.0
    mass = cells["rho"]*vols*\
        snap.info["unit_mass"].express(C.g)
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    etherm = 1.0/(gamma - 1.0) * (mass / mH) * kB * temp
    return np.sum(etherm)

def ekininsnap(snap):
    print "Finding kinetic energy in snap", snap.iout
    amr = snap.amr_source(["rho","vel"])
    cells = CellsToPoints(amr).flatten()
    umass = snap.info["unit_mass"].express(C.g) 
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    ue = umass*uvel**2
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes())**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    ekin = 0.5*rhos*vols*spds**2
    ekin =  np.sum(ekin)*ue
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "KINETIC ENERGY FOUND @ ",time,"Myr:", ekin
    return ekin

def energyinsnap(snap):
    time, etherm = etherminsnap(snap)
    time, ekin = ekininsnap(snap)
    etot = ekin+etherm
    return time, (etherm,ekin,etot)

def maxTinsnap(snap):
    print "Finding max T in snap", snap.iout
    amr = snap.amr_source(["rho","P"])
    cells = CellsToPoints(amr).flatten()
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    maxT = temp.max()
    return maxT

def momentuminsnap(snap,nHlow=0.0,nHhigh=0.0):
    print "Finding momentum in snap", snap.iout
    amr = snap.amr_source(["rho","vel"])
    cells = CellsToPoints(amr).flatten()
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes())**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    pos = cells.points-0.5
    rads = pos+0.0
    dist = np.sqrt(np.sum(pos**2,1))
    for i in range(0,3):
        rads[:,i] /= dist
    rvel = np.sum(rads*vels,1)
    rvelp = rvel+0.0
    rvelp[rvelp < 0.0] = 0.0
    moms = rhos*rvelp*vols
    # Apply threshold?
    if nHlow > 0.0:
        moms = moms[rhos > nHlow]
        rhos = rhos[rhos > nHlow]
    if nHhigh > 0.0:
        moms = moms[rhos < nHhigh]
        rhos = rhos[rhos < nHhigh]
    mom =  np.sum(moms)
    umass = snap.info["unit_mass"].express(C.g)
    mom *= umass*uvel
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "MOMENTUM FOUND @ ",time,"Myr:", mom
    return mom 

def totalmomentuminsnap(snap,nHlow=0.0,nHhigh=0.0):
    print "Finding *total* momentum in snap", snap.iout
    amr = snap.amr_source(["rho","vel"])
    cells = CellsToPoints(amr).flatten()
    uvel = snap.info["unit_velocity"].express(C.cm/C.s)
    rhos = cells["rho"]
    vels = cells["vel"]
    vols = (cells.get_sizes())**3.0
    spds = np.sqrt(np.sum(vels**2.0,1))
    moms = rhos*spds*vols
    # Apply threshold?
    if nHlow > 0.0:
        moms = moms[rhos > nHlow]
        rhos = rhos[rhos > nHlow]
    if nHhigh > 0.0:
        moms = moms[rhos < nHhigh]
        rhos = rhos[rhos < nHhigh]
    mom =  np.sum(moms)
    umass = snap.info["unit_mass"].express(C.g)
    mom *= umass*uvel
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "MOMENTUM FOUND @ ",time,"Myr:", mom
    return mom 

def radiusinsnap(snap):
    print "Finding radius of blast in snap", snap.iout
    amr = snap.amr_source(["rho","P","xHII"])
    cells = CellsToPoints(amr).flatten()
    temp = cells["P"]/cells["rho"]*snap.info["unit_temperature"].express(C.K)
    ion = cells["xHII"]
    pos = cells.points-0.5
    dist = np.sqrt(np.sum(pos**2,1))
    thresh = 0.1 # HACK!
    try:
        maxrad = dist[np.where(ion > thresh)].max()
    except:
        maxrad = 0.0
    maxrad *= snap.info["unit_length"].express(C.pc)
    time = snap.info["time"]*snap.info["unit_time"].express(C.Myr)
    print "RADIUS FOUND @ ",time,"Myr:", maxrad
    return maxrad

def funcinsim(sim,func):
    times = []
    hfunc = Hamu.Algorithm(func)
    Ts = [] # Allows generic formats for non-scalar quantities, etc
    #if nsnsub:
    #    simname = sim.Name()
    #    nsnname = simname.replace("-SN","-NSN")
    #    nsnname = nsnname.replace("-MSN","-NSN")
    #    nsnname = nsnname.replace("-HN","-NSN")
    #    nsnname = nsnname.replace("-NE","-NSN")
    #    nsnname = nsnname.replace("-ME","-NSN")
    #    nsnname = nsnname.replace("-MNE","-NSN")
    #    simnsn = Hamu.Simulation(nsnname)
    for snap in sim.Snapshots():
        time = snaptime.Myr(snap)
        T = hfunc(snap)
        Ts.append(T)
        times.append(time)
    # Get value at tsn
    times = np.array(times)
    # Don't cast Ts in case it's a tuple or whatever
    return times, Ts

def runforfunc(func,simnames):
    name = func.__name__
    plt.clf()
    texts = []
    tlog = False
    ylim = None
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        print "Running for simulation", simname
        #time.sleep(1)
        times, Ts = funcinsim(sim,func)
        if tlog:
            times += 1e-6 # Limit -infinity errors
        col = linestyles.colour(simname)
        line = linestyles.line(simname)
        try:
            maxTs = np.max(Ts)
        except:
            maxTs = 0.0
        if maxTs > 0.0 or not tlog:
            plt.plot(times,Ts,label=sim.Label(),
                     linestyle=line,
                     color=col)
            tx = float(times[-1])
            ty = float(Ts[-1])
            print "TEXT POS", tx, ty
            texts.append(plt.text(tx,ty,simname.replace("_fullrt",""),fontsize="x-small"))
    adjust_text(texts,arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
    # Make Wind property legend
    #if showwindlegend:
    #    plotlines, labels = linestyles.uvlegend()
    #    leg1 = plt.legend(plotlines,labels,loc="upper left",
    #                      ncol=1,fontsize="small",frameon=False)
    #    if showphotonlegend:
    #        plt.gca().add_artist(leg0)
    #    plt.gca().add_artist(leg1)
    #leg0.get_frame().set_edgecolor('w')
    # Quantity-specific stuff
    if tlog:
        plt.xscale("log")
        plt.xlim([1e-6,1])
        
    plt.yscale("log")
    if "maxT" in name:
        plt.ylabel("Temperature / K")
    if "momentum" in name:
        plt.ylabel("Total Momentum / g cm/s")
    if name == "stellarmassinsnap":
        plt.ylabel(r"$M_{*}$ / "+Msolar)
        plt.ylim([10,100])
    # Check for user-specified ylim
    if ylim is not None:
        plt.ylim(ylim)
    if "radius" in name:
        plt.ylabel("radius / pc")
        plt.yscale("linear")
    plt.xlabel("Time / Myr")
    ws = Hamu.CurrentWorkspace()
    suffix = ""#filesuffix
    if suffix == "":
        if tlog:
            suffix += "_tlog"
    plt.savefig("../plots/"+name+"_"+ws+suffix+".pdf")

def prerun(funcs,simnames):
    '''
    Run function once on each snap using only one amr_source call
    Uses deep magic with eval, exec, be warned!
    '''
    funcnames = [f.__name__ for f in funcs]
    # Skim the functions to get the code
    amrs, funccodes = skimfunc.skimfuncs(funcnames,__file__)
    # Overwrite the existing functions
    algs = []
    for name, code in funccodes.iteritems():
        exec(code)
        algs += [Hamu.Algorithm(locals()[name])]
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        for snapold in sim.Snapshots():
            snap = copy.copy(snapold)
            # Make full AMR (bad hack...)
            cached = False
            for alg in algs:
                if not alg.Cached(snap) and not cached:
                    skimfunc.makefullamr(snap.RawData(),amrs)
                    cached = True
                alg(snap)

def run(funcs,simnames):
    '''
    Preloads and runs funcs
    Groups functions for snapshots to minimise memory reloading
    '''
    # Preload the results
    prerun(funcs,simnames)
    # Run the plotter
    for func in funcs:
        runforfunc(func,simnames)

if __name__=="__main__":
    simnames = allsims#[0:12]
    funcs = [stellarmassinsnap,
             maxstellarmassinsnap,
             momentuminsnap,
             totalmomentuminsnap,
             meaniondensinsnap,
             meanionpressinsnap,
             massinsnap]
    run(funcs,simnames)
