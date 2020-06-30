'''
Find the density PDF
Sam Geen, January 2018
'''

from startup import *

from turbustat.statistics import DeltaVariance
from astropy.io import fits
from astropy import units

import columndensity, images

def _densitypdf(snap,los="z"):
    dmap = columndensity.DensityMap(snap,los)
    densities = dmap.Flatten().NH()
    densities = np.log10(densities)
    hist, edges = np.histogram(densities,bins=20)
    centres = 0.5*(edges[1:] + edges[:-1])
    hist = hist/float(densities.shape[0])
    return centres, hist
densitypdf = Hamu.Algorithm(_densitypdf)

def runforsim(simname):
    print "Running for simulation", simname     
    sim = hamusims[simname]
    
    tcodes = tffcloud_code * np.array([0.5,1.0])
    lines = ["--","-"]
    labels = ["0.5 $t_{ff}$","$t_{ff}$"]
    plt.clf()
    for line, tcode, label in zip(lines, tcodes,labels):
        snap = sim.FindAtTime(tcode)
        centres, hist = densitypdf(snap,los="z")
        plt.plot(centres, hist, linestyle=line,label=label)
    plt.xlabel("log($N_{H}$ / cm$^{-2}$)")
    plt.ylabel("Density PDF")
    plt.yscale("log")
    plt.legend(frameon=False,fontsize="x-small")
    plt.savefig("../plots/pdfdensity.pdf")

def run(simnames,plotname):
    plt.clf()
    for simname in simnames:
        t, slope = runforsim(simname)
        import pdb; pdb.set_trace()

if __name__=="__main__":
    run(imfsims,"imf")
    #run(icsims,"ic")
