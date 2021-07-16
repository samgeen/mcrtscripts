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
    mask = np.where(hydro.T[0:nx] > 1e3)[0]
    if len(mask) > 0:
        radius = hydro.x[mask[-1]]
    return radius
hydrofuncs["shellradius"] = findshellradius
hydrolabels["shellradius"] = "Shell Radius / pc"
hydrolabelsdiff["shellradius"] = "Shell Velocity / km/s"
hydrounits["shellradius"] = wunits.pc
hydrounitsdiff["shellradius"] = 1e5 # km/s in cgs

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

def plottimefuncs(sims,hydro,diff=False,yscale="log"):
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
    if hydro == "shellradius" and diff:
        plt.plot(time/tunit,41.4 * (0.2 / (5.0/11.0)**0.33))
    plt.xlabel("Time / kyr")
    plt.ylabel(label)
    plt.legend(fontsize="x-small",frameon=False)
    plt.xscale("log")
    plt.yscale(yscale)
    plt.savefig("../plots/oned/timefunc_"+hydro+difftxt+".pdf")

def plotprofiles(snap):
    plt.clf()
    simname = snap.Name()
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    r = hydro.x[0:nx] / wunits.pc
    vol = hydro.vol[0:nx]
    # Density
    plt.plot(r,hydro.nH[0:nx],label="$n_\mathrm{H}$ / cm$^{-3}$")
    # Temperature
    plt.plot(r,hydro.T[0:nx],label="T / K")
    # Cooling
    dTdt = weltgeist.cooling.TemperatureChange(1.0)
    dEdt = 1.5 * hydro.nH[0:nx] * wunits.kB * dTdt
    plt.plot(r,-dEdt[0:nx]*1e24,label="Cooling x1e24 / erg/s/cm$^3$")
    # Entropy
    entropy = findentropy(snap,donotsum=True)
    plt.plot(r,entropy,label="Entropy / K cm$^{2}$")
    # Finish plotting
    plt.legend(fontsize="x-small",frameon=False)
    plt.xlabel("Radius / pc")
    plt.yscale("log")
    os.makedirs("../plots/oned/"+simname, exist_ok=True)
    plt.savefig("../plots/oned/"+simname+"/1Dhydrotest_"+str(snap.OutputNumber())+".pdf")
        

def CompareSims(sims):
    """
    Compare two sets of sims
    """
    # Plot hydro variables over time
    print("---")
    print("TIME FUNCS")
    print("---")
    for hydro in hydrofuncs.keys():
        print (hydro)
        plottimefuncs(sims,hydro)
    plottimefuncs(sims,"shellradius",diff=True,yscale="linear")
    # Plot profiles
    print("---")
    print("PROFILES")
    print("---")
    for simname, sim in sims.items():
        print(simname)
        reload(weltgeist)
        for snap in sim.Snapshots():
            print (snap.OutputNumber())
            plotprofiles(snap)


if __name__=="__main__":
    filenames = ["../outputs/WindOnly1e+36ergspersecond1000_w2_N2048"]
    simnames = ["Wind-L36"]
    sims = {}
    for simname, filename in zip(simnames, filenames):
        sims[simname] = Hamu.MakeSimulation(simname,filename)
    CompareSims(sims)

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