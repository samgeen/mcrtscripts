'''
Find the total star formation efficiency
Sam Geen, December 2017
'''

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

from startup import *

import matplotlib.patheffects as pe

def _tsfeinsnap(snap):
    sink = sinks.FindSinks(snap)
    mgas = 1e4
    return np.sum(sink.mass) / mgas

tsfeinsnap = Hamu.Algorithm(_tsfeinsnap)

def runforsim(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    return timefuncs.timefunc(sim,tsfeinsnap)

def run(simnames,plotname):
    plt.clf()
    for simname in simnames:
        t, sfe = runforsim(simname)
        plt.plot(t,sfe,color=linestyles.Colour(simname),label=linestyles.Label(simname),
                 path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
    plt.xlabel("Time / Myr")
    plt.ylabel("SFE $\equiv M_{stars} / M_{ini}$")
    plt.legend(frameon=False,fontsize="x-small")
    plt.ylim([-0.01,0.25])
    plt.savefig(plotfolder+"tsfe"+plotname+".pdf")
        
if __name__=="__main__":
    run(imfsims,"imf")
    run(icsims, "ic")
