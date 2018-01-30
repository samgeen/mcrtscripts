'''
Calculate escape fraction from ray profiles
Sam Geen, April 2015
'''

import customplot
import matplotlib.pyplot as plt

import numpy as np

import pymses
from pymses.utils import constants as C

import Hamu, rayprof

nray = rayprof.nray
nprofs = rayprof.nprofs

# ionlim = some value where we assume everything is ionised
def FEsc(snap,ionlim=0.5):
    r, p = rayprof.makeprofsHamu(snap.hamusnap,"xHII")
    radii = r[0:nray]
    profs = np.reshape(p,(nprofs, nray)).T
    nesc = 0.0
    for i in range(0,nprofs):
        if np.min(profs[:,i]) > ionlim:
            nesc += 1.0
    fesc = nesc / nprofs
    return fesc

FEscHamu = Hamu.Algorithm(FEsc)

def RunForSim(simname,tstart):
    sim = Hamu.Simulation(simname)
    fescs = []
    myr = None
    for snap in sim.Snapshots():
        if myr is None:
            myr = snap.RawData().info["unit_time"].express(C.Myr)
        fescs.append(FEscHamu(snap))
    times = sim.Times()*myr
    times -= tstart
    plt.plot(times,fescs,label=sim.Label()+" t$_{\mathrm{ff}}$="+str(tstart)+" Myr")
    

def Run(simnames,tstarts):
    plt.clf()
    for simname in simnames:
        RunForSim(simname,tstarts[simname])
    plt.xlim([0,4])
    plt.legend(fontsize=8,loc="lower left")
    plt.xlabel("Time / Myr")
    plt.ylabel("Escape Fraction")
    plt.yscale("log")
    plt.ylim([0.01,1])
    plt.savefig("../plots/escapefractions/escapefractions_all.pdf")

if __name__=="__main__":
    tstarts = {}
    tstarts["N47_M4_B02"] = 1.25
    tstarts["N48_M4_B02"] = 1.25
    tstarts["N49_M4_B02"] = 1.25
    tstarts["N48_M4_B02_C"] = 1.25*0.5**3
    tstarts["N48_M4_B02_C2"] = 1.25*0.75**3
    simnames = ["N47_M4_B02","N48_M4_B02","N49_M4_B02",
                "N48_M4_B02_C","N48_M4_B02_C2"]
    Run(simnames,tstarts)
    
