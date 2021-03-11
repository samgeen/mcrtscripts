'''
Find parameters in a hydro snapshot in Weltgeist
Sam Geen, January 2021
'''

import os, sys

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

hydrofuncs = {}
hydrolabels = {}
hydrounits = {}
def findshellradius(snap):
    """
    Find the radius of the shell
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
hydrounits["shellradius"] = wunits.pc

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

def findcooling(snap,maskon=True):
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    dTdt = weltgeist.cooling.TemperatureChange(1.0)
    if maskon and "maskCD" in snap.Folder():
        dTdt = weltgeist.cooling.MaskContactDiscontinuity(dTdt)
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
        dTdt = weltgeist.cooling.MaskContactDiscontinuity(dTdt)
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
    dTdt2 = weltgeist.cooling.MaskContactDiscontinuity(dTdt)
    # Recover the CD part
    CDcooling = dTdt - dTdt2
    import pdb; pdb.set_trace()
    dEdt = 1.5 * hydro.nH[0:nx] * hydro.vol[0:nx] * wunits.kB * CDcooling
    # Remove heating
    dEdt[dEdt < 0.0] = 0.0
    return np.sum(dEdt)
hydrofuncs["CDcooling"] = findCDcooling
hydrolabels["CDcooling"] = "Cooling rate in the CD / erg / s"
hydrounits["CDcooling"] = 1.0 # cgs already

def plottimefuncs(sims,hydro):
    plt.clf()
    tunit = wunits.year/1e3
    yunit = hydrounits[hydro]
    for simname, sim in sims.items():
        time, yvals = timefuncs.timefunc(sim,hydrofuncs[hydro])
        plt.plot(time/tunit,yvals/yunit,label=simname)
    plt.xlabel("Time / kyr")
    plt.ylabel(hydrolabels[hydro])
    plt.legend(fontsize="x-small",frameon=False)
    plt.yscale("log")
    plt.savefig("../plots/oned/timefunc_"+hydro+".pdf")

def plotprofiles(snap):
    plt.clf()
    simname = snap.Name()
    hydro = snap.RawData().hydro
    nx = hydro.ncells
    r = hydro.x[0:nx] / wunits.pc
    plt.plot(r,hydro.nH[0:nx],label="$n_\mathrm{H}$ / cm$^{-3}$")
    plt.plot(r,hydro.T[0:nx],label="T / K")
    dTdt = weltgeist.cooling.TemperatureChange(1.0)
    dEdt = 1.5 * hydro.nH[0:nx] * wunits.kB * dTdt
    plt.plot(r,-dEdt[0:nx]*1e24,label="Cooling x1e24 / erg/s/cm$^3$")
    plt.legend(fontsize="x-small",frameon=False)
    plt.xlabel("Radius / pc")
    plt.yscale("log")
    try:
        os.mkdir("../plots/")
    except:
        pass
    try:
        os.mkdir("../plots/oned/")
    except:
        pass
    try:
        os.mkdir("../plots/oned/"+simname)
    except:
        pass
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
    # Plot profiles
    print("---")
    print("PROFILES")
    print("---")
    for simname, sim in sims.items():
        print(simname)
        for snap in sim.Snapshots():
            print (snap.OutputNumber())
            plotprofiles(snap)



filenames = ["../outputs/35Msun_n1000_w2_N256_uvwind_coolingfix2maskCD/","../outputs/35Msun_n1000_w2_N256_uvwind_coolingfix2/"][::-1]
simnames = ["No CD Cooling","CD Cooling"][::-1]
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