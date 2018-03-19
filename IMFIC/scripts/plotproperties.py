'''
Find the total star formation efficiency
Sam Geen, December 2017
'''

from pymses.filters import CellsToPoints
from pymses.utils import constants as C

from startup import *

import matplotlib.patheffects as pe

import findproperties, starrelations

def _tsfeinsnap(snap):
    sink = sinks.FindSinks(snap)
    mgas = 1e4
    return np.sum(sink.mass) / mgas

tsfeinsnap = Hamu.Algorithm(_tsfeinsnap)

def tsfe(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    return timefuncs.timefunc(sim,tsfeinsnap)

momentuminsnap = Hamu.Algorithm(findproperties.totalmomentuminsnap)
def momentum(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    return timefuncs.timefunc(sim,momentuminsnap)

radiusinsnap = Hamu.Algorithm(findproperties.radiusinsnap)
def radius(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    return timefuncs.timefunc(sim,radiusinsnap)

def _NphotonsHII(snap):
    return np.sum(starrelations.findNphotons(snap.hamusnap)[3:5])
NphotonsHII = Hamu.Algorithm(_NphotonsHII)

def nphotonsHII(simname):
    print "Running for simulation", simname
    sim = hamusims[simname]
    return timefuncs.timefunc(sim,NphotonsHII)

def run(simfunc,simnamesets,plotlabels):
    plt.clf()
    funcname = simfunc.__name__
    fig, axes = plt.subplots(1,2,sharex=True,sharey=True)
    first = True
    for ax, simnames, plotlabel in zip(axes,simnamesets,plotlabels):
        for simname in simnames:
            t, y = simfunc(simname)
            ax.plot(t,y,color=linestyles.Colour(simname),label=linestyles.Label(simname),
                     path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
        ax.set_xlabel("Time / Myr")
        if first:
            if funcname == "tsfe":
                ax.set_ylabel("SFE $\equiv M_{stars} / M_{ini}$")
            if funcname == "momentum":
                ax.set_ylabel("Momentum / g cm / s")
            if funcname == "radius":
                ax.set_ylabel("Mean HII Region radius / pc")
            if funcname == "nphotonsHII":
                ax.set_ylabel("Ionising Photon Emission Rate / s$^{-1}$")
        ax.legend(frameon=False,fontsize="x-small")
        #ax.set_ylim([-0.01,0.25])
        ax.text(0.95, 0.95,plotlabel, ha='right', va='top', transform=ax.transAxes)
        first = False
    fig.subplots_adjust(wspace=0)
    fig.set_size_inches(14,6)
    fig.savefig(plotfolder+funcname+"_both.pdf", dpi=80)
        
if __name__=="__main__":
    for func in [nphotonsHII,tsfe,momentum,radius]:
        run(func,(imfsims,icsims),(linestyles.starsc,linestyles.turbsc))
