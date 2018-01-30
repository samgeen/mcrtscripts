'''
Find the distribution of radii to the diffuse medium from the source
Sam Geen, February 2015
'''

import customplot, rayprof, profilesphere
import pymses, Hamu
import os
import numpy as np
import matplotlib.pyplot as plt
from pymses.utils import constants as C

rextent = 0.5 # fraction of box *radius* (1/2 boxlen)
pctocm = 3.08567758e18
hcut = 1e4 # density cut
def FindFitExp(r,p,rcut=1e10,params=True):
    p = p[r < rcut]
    r = r[r < rcut]
    r = r[1:]
    p = p[1:]
    notnan = np.logical_not(np.isnan(r*p))
    r = r[notnan]
    p = p[notnan]
    p[p == 0.0] = p[p > 0.0].min() # Safety net
    lnp = np.log(p)
    lnr = np.log(r)
    fit = np.polyfit(r,lnp,deg=1)
    #     rho = rho0 * exp(-r/r0)
    # --> ln(rho) = -r/r0 + ln(rho0)
    mult,const = fit
    rho0 = np.exp(const)
    r0 = -1.0/mult
    print "rho0, r0:", rho0, r0
    pfit = rho0 * np.exp(-r/r0)
    if params:
        return rho0, r0
    else:
        return r,pfit

def FindFitPow(r,p,rcut=1e10,params=True):
    p = p[r < rcut]
    r = r[r < rcut]
    r = r[1:]
    p = p[1:]
    notnan = np.logical_not(np.isnan(r*p))
    r = r[notnan]
    p = p[notnan]
    p[p == 0.0] = p[p > 0.0].min() # Safety net
    lnp = np.log(p)
    lnr = np.log(r)
    fit = np.polyfit(lnr,lnp,deg=1)
    #     rho = rho0 * r ^ (-n)
    # --> ln(rho) = -n * ln(r) + ln(rho0)
    mult,const = fit
    rho0 = np.exp(const)
    n = -mult
    print "rho0, n:", rho0, n
    pfit = rho0 * r ** (-n)
    if params:
        return rho0, n
    else:
        return r,pfit


FindFit = FindFitExp


def Plot(simname,time):
    print "Running testcutprofile for ", simname
    # Find snapshot
    sim = Hamu.Simulation(simname)
    first = sim.Snapshots()[0]
    tcode = time / first.RawData().info["unit_time"].express(C.Myr)
    snap = sim.FindAtTime(tcode)
    # Run
    boxrad = snap.RawData().info["boxlen"]/2.0
    r,p = rayprof.cutprofileHamu(snap,"rho",hcut) 
    # Plot
    plt.clf()
    plt.plot(r,p)
    if np.max(p) > 0.0:
        rf,pf = FindFit(r,p,rcut=boxrad*rextent,params=False)
        plt.plot(rf,pf,"--")
    plt.yscale("log")
    #plt.xscale("log")
    plt.savefig("../plots/testcutprofile/testcutprofile_"+simname+".pdf")
    print "Done for ", simname
    print "-----------------"

if __name__=="__main__":
    Plot("N48_M4_B02",1.25)
    Plot("N48_M4_B02_F2",1.25*2)
    Plot("N48_M4_B02_F3",1.25*3)
    Plot("N48_M4_B02_C",1.25*0.5**3)
    Plot("N48_M4_B02_C2",1.25*0.75**3)