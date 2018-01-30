'''
Find the PDF of maximum density along a line of sight
Sam Geen, May 2015
'''

import sys
sys.path.append("/home/sgeen/MC_RT/scripts")

import customplot
import matplotlib.pyplot as plt

from pymses.utils import constants as C

import Hamu
import rayprof, rdotplot

import linestyles

import numpy as np

def findmaxes(snap):
    r, profs = rayprof.makeprofsHamu(snap.hamusnap,"rho")
    ninray = rayprof.nray
    nprofs = rayprof.nprofs
    profs = np.reshape(profs,(nprofs,ninray))
    # shape of profs is (ninray, nprofs)
    maxes = np.amax(profs,axis=1)
    maxes = np.sort(maxes)
    nums = np.arange(1,nprofs+1)
    return nums, maxes

findmaxesHamu = Hamu.Algorithm(findmaxes)

def runforsims(simnames,labels,tstart):
    plt.clf()
    for simname, label in zip(simnames,labels):
        sim = Hamu.Simulation(simname)
        utime = sim.Snapshots()[0].RawData().info["unit_time"].express(C.Myr)
        start = sim.FindAtTime(tstart/utime)
        nums, maxes = findmaxesHamu(start)
        nums = nums/float(rayprof.nprofs)
        #nums = np.log10(nums)
        plt.plot(maxes[::-1],nums,
                 color=linestyles.col(simname),
                 label=label)
    plt.xlabel("Maximum n$_{\mathrm{H}}$ along line of sight / cm$^{-3}$")
    plt.ylabel("P($>$n$_{\mathrm{H}}$)")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize="small",frameon=False)
    plt.savefig("../plots/losmax/losmax.pdf")

if __name__=="__main__":
    tsn = 5.52990445
    simnames = ["N00-NSN","N49-NSN","N50-NSN","N51-NSN"]
    labels = [linestyles.emlabel(sim) for sim in simnames]
    runforsims(simnames,labels,tstart=tsn)
