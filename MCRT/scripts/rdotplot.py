'''
Plot rdot against time for sims
Sam Geen, May 2015
'''

import customplot

import matplotlib.pyplot as plt
import numpy as np

import Hamu
import outflowmodel

def Findtstart(simname):
    tstart = 1.25
    if "_C2" in simname:
        tstart *= 0.75**3
    elif "_C" in simname:
        tstart *= 0.5**3
    return tstart

def run(simnames,name=""):
    plt.clf()
    cs = outflowmodel.cs/1e5
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        times, rdots = outflowmodel.FindRdotPowerLaw(sim)
        plt.plot(times,rdots,label=sim.Label())
        # TODO: DIFFERENT LINE COLOURS, STYLES
    plt.plot(times,np.zeros((len(times)))+cs)
    plt.xlim([0,5])
    plt.xlabel("Times / Myr")
    plt.ylabel("Shell Velocity / km/s")
    plt.legend(fontsize="xx-small")
    if len(name) > 0:
        name = "_"+name
    plt.savefig("../plots/rdot/rdot"+name+".pdf")

# Trying something else here
def roverrs(simnames,name=""):
    plt.clf()
    for simname in simnames:
        tstart = Findtstart(simname)
        outflowmodel.tstart = tstart
        sim = Hamu.Simulation(simname)
        times, riis = outflowmodel.FindriiPowerLaw(sim)
        rs = outflowmodel.Findrstromgren(sim)/outflowmodel.pcincm
        index = outflowmodel.FindProfilePower(sim)
        powfact = 3.0/4.0 + index/2.0
        riis /= rs
        riis = riis**powfact
        times /= tstart
        plt.plot(times,riis,label=sim.Label())
        print "Power law index:", index
        # TODO: DIFFERENT LINE COLOURS, STYLES
    plt.xlim([0,5])
    plt.yscale("log")
    plt.xlabel("Times / $t_{ff}$")
    plt.ylabel("(Radius / R$_S$)^(3/4-w/2)")
    plt.legend(fontsize="xx-small")
    if len(name) > 0:
        name = "_"+name
    plt.savefig("../plots/rdot/riis"+name+".pdf")

if __name__=="__main__":
    simnames = ["N47_M4_B02","N48_M4_B02","N49_M4_B02",
                "N48_M4_B02_C","N48_M4_B02_C2"]
    roverrs(simnames,"all")
