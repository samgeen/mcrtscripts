'''
Compare energy in runs
Sam Geen, April 2015
'''

import customplot
import matplotlib.pyplot as plt

import pymses
from pymses.utils import constants as C

import numpy as np

import Hamu
import profilesphere
profilemodule = profilesphere

Te = 8400.0
kB = 1.3806e-16 # erg / K
gamma = 1.4
X = 0.74
#mu = X*2.0 + 0.25*(1-X)*2.0 # Ionised hydrogen plus once-ionised He
mu = 0.61 # From Matzner 2002
mH = 1.67e-24 # g
cs = np.sqrt(gamma * kB * Te / (mH*mu))
print "USING cs = ", cs/1e5, "km/s"
G = 6.674e-8
pcincm = 3.08567758e18
Myrins = 3.15569e13
kmincm = 1e5

def plotenergies(simname):
    sim = Hamu.Simulation(simname)
    KEs = []
    GPEs = []
    times = []
    myr = None
    for snap in sim.Snapshots():
        if myr is None:
            myr = snap.RawData().info["unit_time"].express(C.Myr)
        # Get profiles and scale to cgs
        r,rho = profilemodule.profileHamu(snap,"rho",1e6)
        r,spd = profilemodule.profileHamu(snap,"spd",1e6)
        r *= pcincm
        rho *= mH / X
        spd *= kmincm
        # ASSUME THAT r IS UNIFORMLY SPACED!!!
        dr = r[1] - r[0]
        # r SHOULD BE CENTRED ON THE SAMPLED SPHERICAL SHELL
        mass = rho * r**2 * dr
        cmass = np.cumsum(mass)
        KE = 0.5 * spd**2 * mass
        GPE = G * cmass * mass / r
        KEs.append(np.sum(KE))
        GPEs.append(np.sum(GPE))
        times.append(snap.Time()*myr)
    times = np.array(times)
    GPEs = np.array(GPEs)
    KEs = np.array(KEs)
    ratios = KEs / GPEs
    plt.plot(times,KEs,label="KE")
    plt.plot(times,GPEs,label="GPE")
    plt.xlabel("Time / Myr")
    plt.ylabel("Energy / erg")
    plt.yscale("log")
    plt.legend()
    ax2 = plt.gca().twinx()
    ax2.plot(times,ratios,"k")
    ax2.set_ylabel("KE / GPE")
    ax2.set_yscale("log")
    plt.savefig("../plots/energycompare/energycompare_"+simname+".pdf")

if __name__=="__main__":
    simnames = ["N00_M4_B02"]
    for simname in simnames:
        plotenergies(simname)
