'''
Find parameters in a hydro snapshot in Weltgeist
Sam Geen, January 2021
'''

import os, sys
from importlib import reload  

import numpy as np

sys.path.append("../../scripts") 
import HamuLite as Hamu
import timefuncs

sys.path.append("/home/samgeen/Programming/Astro/Weltgeist")
import weltgeist
import weltgeist.units as wunits

sys.path.append("/home/samgeen/Programming/Astro/WindInUV")
import stars

import matplotlib.pyplot as plt

import snaptime
snaptime.nopymses = True

hydrofuncs = {}
hydrolabels = {}
hydrolabelsdiff = {}
hydrounits = {}
hydrounitsdiff = {}
def findshellradius(snap):
    """
    Find the inner radius of the shell
    """
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    radius = 0.0
    mask = np.where(hydro.T[0:nx] > 5e4)[0]
    if len(mask) > 0:
        radius = hydro.x[mask[-1]]
    return radius
hydrofuncs["shellradius"] = findshellradius
hydrolabels["shellradius"] = "Shell Radius / pc"
hydrolabelsdiff["shellradius"] = "Shell Velocity / km/s"
hydrounits["shellradius"] = wunits.pc
hydrounitsdiff["shellradius"] = 1e5 # km/s in cgs


def findhiiradius(snap):
    """
    Find the HII region radius
    """
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    radius = 0.0
    mask = np.where(hydro.xhii[0:nx] > 0.1)[0]
    if len(mask) > 0:
        radius = hydro.x[mask[-1]]
    return radius
hydrofuncs["hiiradius"] = findhiiradius
hydrolabels["hiiradius"] = "Ionisation Front Radius / pc"
hydrolabelsdiff["hiiradius"] = "Ionisation Front Velocity / km/s"
hydrounits["hiiradius"] = wunits.pc
hydrounitsdiff["hiiradius"] = 1e5 # km/s in cgs

def findenergy(snap):
    """
    Find the kinetic energy in the simulation
    """
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    kinetic = np.sum(hydro.KE[0:nx])
    thermal = np.sum(hydro.TE[0:nx])
    return kinetic + thermal
hydrofuncs["energy"] = findenergy
hydrolabels["energy"] = "Total Energy / erg"
hydrounits["energy"] = 1.0 # cgs already

def findhotgasenergy3(snap):
    """
    Find the kinetic energy in the simulation
    """
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    kinetic = hydro.KE[0:nx]
    thermal = hydro.TE[0:nx]
    energy = kinetic + thermal
    mask = hydro.T[0:nx] > 1e5
    energy = np.sum(energy[mask])
    return energy
hydrofuncs["hotgasenergy"] = findhotgasenergy3
hydrolabels["hotgasenergy"] = "Energy($T_{\mathrm{gas}} > 10^5~$K) / erg"
hydrounits["hotgasenergy"] = 1.0 # cgs already

def findentropy(snap,donotsum=False):
    """
    Find the kinetic energy in the simulation
    """
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    temperature = hydro.T[0:nx]
    nH = hydro.nH[0:nx]
    #entropy = 3.0 * wunits.kB * temperature * (rho)**(-2.0/3.0) * (wunits.X / wunits.mH)
    entropy = temperature * (nH)**(-2.0/3.0)
    if not donotsum:
        entropy = np.sum(entropy)
    return entropy
hydrofuncs["entropy"] = findentropy
hydrolabels["entropy"] = "Total Entropy / K cm$^{2}$"
hydrounits["entropy"] = 1.0 # cgs already

def findcooling(snap,maskon=True):
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    dTdt = weltgeist.cooling.TemperatureChange(1.0)
    if maskon and "maskCD" in snap.Folder():
        if "v2" in snap.Folder():
            dTdt = weltgeist.cooling.MaskContactDiscontinuityV2(dTdt)
        else:
            dTdt = weltgeist.cooling.MaskContactDiscontinuityV1(dTdt)
    # Find Cooling (so -ve temperature change)
    dEdt = - 1.5 * hydro.nH[0:nx] * hydro.vol[0:nx] * wunits.kB * dTdt
    # Remove heating
    dEdt[dEdt < 0.0] = 0.0
    return np.sum(dEdt)
hydrofuncs["cooling"] = findcooling
hydrolabels["cooling"] = "Cooling rate / erg / s"
hydrounits["cooling"] = 1.0 # cgs already

def findhotgascooling(snap,maskon=True):
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    dTdt = weltgeist.cooling.TemperatureChange(1.0)
    if maskon and "maskCD" in snap.Folder():
        if "v2" in snap.Folder():
            dTdt = weltgeist.cooling.MaskContactDiscontinuityV2(dTdt)
        else:
            dTdt = weltgeist.cooling.MaskContactDiscontinuityV1(dTdt)
    # Find Cooling (so -ve temperature change)
    dEdt = - 1.5 * hydro.nH[0:nx] * hydro.vol[0:nx] * wunits.kB * dTdt
    mask = hydro.T[0:nx] < 1e5
    dEdt[mask] = 0.0
    # Remove heating
    dEdt[dEdt < 0.0] = 0.0
    return np.sum(dEdt)
hydrofuncs["hotgascooling"] = findhotgascooling
hydrolabels["hotgascooling"] = "Cooling rate gas above $10^5$ K / erg / s"
hydrounits["hotgascooling"] = 1.0 # cgs already

def findCDcooling(snap,maskon=True):
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    # Find Cooling (so -ve temperature change)
    dTdt = -weltgeist.cooling.TemperatureChange(1.0)
    # Find cooling minus the CD part
    dTdt2 = weltgeist.cooling.MaskContactDiscontinuityV2(dTdt)
    # Recover the CD part
    CDcooling = dTdt - dTdt2
    dEdt = 1.5 * hydro.nH[0:nx] * hydro.vol[0:nx] * wunits.kB * CDcooling
    # Remove heating
    dEdt[dEdt < 0.0] = 0.0
    return np.sum(dEdt)
hydrofuncs["CDcooling"] = findCDcooling
hydrolabels["CDcooling"] = "Cooling rate in the CD / erg / s"
hydrounits["CDcooling"] = 1.0 # cgs already

def plottimefuncs(setname,sims,hydro,diff=False,xscale="log",yscale="log"):
    plt.clf()
    tunit = wunits.year*1e3
    yunit = hydrounits[hydro]
    label = hydrolabels[hydro]
    difftxt = ""
    for simname, sim in sims.items():
        reload(weltgeist)
        time, yvals = timefuncs.timefunc(sim,hydrofuncs[hydro])
        if diff:
            yvals = np.diff(yvals) / np.diff(time)
            time = 0.5*(time[1:]+time[:-1])
            yunit = hydrounitsdiff[hydro]
            label = hydrolabelsdiff[hydro]
            difftxt = "_diff"
        plt.plot(time/tunit,yvals/yunit,label=simname)
        #print(yvals.min(), yvals.max(),hydro)
    # if hydro == "shellradius" and diff:
    #     if len(time) > 0:
    #         plt.plot(time/tunit,41.4 * (0.2 / (5.0/11.0)**0.33)+time*0.0)
    plt.xlabel("Time / kyr")
    plt.ylabel(label)
    plt.legend(fontsize="x-small",frameon=False)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.savefig("../plots/oned/timefunc_"+setname+"_"+hydro+difftxt+".pdf")

def plotprofiles(setname,snap):
    plt.clf()
    simname = snap.Name()
    hydro = snap.RawData().hydro
    time = snap.Time() / wunits.year / 1e3
    nx = hydro.ncells
    r = hydro.x[0:nx] / wunits.pc
    vol = hydro.vol[0:nx]
    # Density
    plt.plot(r,hydro.nH[0:nx],label="$n_\mathrm{H}$ / cm$^{-3}$")
    # Temperature
    plt.plot(r,hydro.T[0:nx],label="T / K")
    # Ionisation fraction
    plt.plot(r,hydro.xhii[0:nx],label="$x_{\mathrm{HII}}$")
    # Cooling
    dTdt = weltgeist.cooling.TemperatureChange(1.0)
    dEdt = 1.5 * hydro.nH[0:nx] * wunits.kB * dTdt
    plt.plot(r,-dEdt[0:nx]*1e24,label="Cooling x1e24 / erg/s/cm$^3$")
    # Entropy
    #entropy = findentropy(snap,donotsum=True)
    #plt.plot(r,entropy,label="Entropy / K cm$^{2}$")
    # Finish plotting
    plt.legend(fontsize="x-small",frameon=False,title=("%.2f" % time)+" kyr")
    plt.xlabel("Radius / pc")
    plt.yscale("log")
    os.makedirs("../plots/oned/"+setname, exist_ok=True)
    os.makedirs("../plots/oned/"+setname+"/"+simname, exist_ok=True)
    plt.savefig("../plots/oned/"+setname+"/"+simname+"/1Dhydrotest_"+str(snap.OutputNumber()).zfill(5)+".pdf")

    
def plotenergyretained(setname,sims,starmass,starmetal,starrotating=True):
    plt.clf()
    tunit = wunits.year*1e3
    difftxt = ""
    star = stars.Star(starmass,starmetal,starrotating)
    for simname, sim in sims.items():
        reload(weltgeist)
        times, gasenergy = timefuncs.timefunc(sim,findhotgasenergy3)
        windenergies = np.array([star.WindEnergy(t) for t in times])
        yvals = gasenergy / windenergies
        plt.plot(times/tunit,yvals,label=simname)
    plt.xlabel("Time / kyr")
    plt.ylabel("$E_{hot} / E_{wind}$")
    plt.legend(fontsize="x-small",frameon=False)
    plt.plot([0.1,200],[1.0/3.0]*2,"k--",alpha=0.5)
    plt.xlim([0.1,200])
    plt.xscale("linear")
    plt.yscale("linear")
    plt.savefig("../plots/oned/timefunc_"+setname+"_energyretained"+difftxt+".pdf")
        

def CompareSims(setname,sims):
    """
    Compare sets of simulations
    """
    # Plot hydro variables over time
    print("---")
    print("TIME FUNCS")
    print("---")
    for hydro in hydrofuncs.keys():
        print (hydro)
        plottimefuncs(setname,sims,hydro)
    plottimefuncs(setname,sims,"shellradius",diff=True,xscale="linear",yscale="log")
    plottimefuncs(setname,sims,"hiiradius",diff=True,xscale="linear",yscale="log")
    # Plot profiles
    print("---")
    print("PROFILES")
    print("---")
    for simname, sim in sims.items():
        print(simname)
        reload(weltgeist)
        for snap in sim.Snapshots():
            print (snap.OutputNumber())
            plotprofiles(setname,snap)



def MakeSims(simnames,filenames):
    sims = {}
    for simname, filename in zip(simnames, filenames):
        sims[simname] = Hamu.MakeSimulation(simname,filename,simname)
    return sims

if __name__=="__main__":
    # Step function test
    simnames = ["Step Test"]
    filenames = ["../outputs/35Msun_n100_w2_N2048_uvwind_coolingfix2_Bfield_stepout"]
    stepset = MakeSims(simnames,filenames)
    # Feedback physics comparison
    simnames = ["UV only","UV+wind","UV+wind+maskCD"]
    filenames = ["../outputs/35Msun_n100_w2_N2048_uvonly_coolingfix2maskCD_v2_Bfield",
                 "../outputs/35Msun_n100_w2_N2048_uvwind_coolingfix2_Bfield",
                 "../outputs/35Msun_n100_w2_N2048_uvwind_coolingfix2maskCD_v2_Bfield"]
    fbset = MakeSims(simnames,filenames)
    # CD Mask comparison
    simnames = ["No cooling", "CD Mask","No CD Mask"]
    filenames = ["../outputs/35Msun_n100_w2_N2048_windonly_coolingfix2_nocooling",
                 "../outputs/35Msun_n100_w2_N2048_windonly_coolingfix2maskCD_v2_Bfield",
                 "../outputs/35Msun_n100_w2_N2048_windonly_coolingfix2_Bfield"]
    windonlyset = MakeSims(simnames,filenames)
    # Density comparison
    simnames = ["$10^3~$cm$^{-3}$","$10^2~$cm$^{-3}$"]
    filenames = ["../outputs/35Msun_n100_w2_N2048_uvonly_coolingfix2maskCD_v2_Bfield",
                 "../outputs/35Msun_n100_w2_N2048_uvwind_coolingfix2_Bfield"]
    densityset = MakeSims(simnames,filenames)
    # Big list of simulation sets
    simsets = {"stepset":stepset,"fbset":fbset,"windonly_cdmaskset":windonlyset,"densityset":densityset}

    # Plot energy retained
    plotenergyretained("windonly_cdmaskset",windonlyset,35.0,0.014,True)

    # Make other general plots
    for setname, sims in simsets.items():
        CompareSims(setname,sims)

'''
#simname = "../outputs/30Msun_n3000_w2_N128_windonly_coolingfix/"
#simname = "../outputs/30Msun_n3000_w2_N256_uvwind_coolingfix2"
#simname = "../outputs/35Msun_n1000_w2_N256_uvwind_coolingfix2"
# CONTACT DISCONTINUITY COOLING TURNED OFF
#simname = "../outputs/35Msun_n1000_w2_N256_uvwind_coolingfix2maskCD/"
# CD COOLING TURNED ON
simname = "../outputs/35Msun_n1000_w2_N256_uvwind_coolingfix2/"
#sim = Hamu.MakeSimulation("TEST","../outputs/30Msun_n3000_w2_N512_uvwind/")
#sim = Hamu.MakeSimulation("TEST","../outputs/test/")
sim = Hamu.MakeSimulation("TEST",simname)
times, radii = timefuncs.timefunc(sim,findshellradius)
vels = np.diff(radii)/np.diff(times)
avtimes = 0.5*(times[1:] + times[:-1])
print(times,radii)
print(avtimes,vels/1e5)
times, energies = timefuncs.timefunc(sim,findenergy)
times, cooling = timefuncs.timefunc(sim,findcooling)
#print(times,cooling)
star = stars.Star(30,0.014,rotating=True)
windenergy = star.WindEnergy(times[-1])
print(times,energies)
print(windenergy)
print(energies[-1] / windenergy)
for snap in sim.Snapshots()[::-1]:
    plotcooling(snap)
'''