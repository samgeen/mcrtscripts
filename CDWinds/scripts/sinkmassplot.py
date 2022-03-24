'''
Plot the main star's sink mass as a function of time
Sam Geen, March 2022
'''

from startup import *

import matplotlib.patheffects as pe

import sinks, stellars

def _findmainsinkmass(snap):
    stellar = stellars.FindStellar(snap)
    if len(stellar.mass) == 0:
        print("No stellar objects, returning 0.0")
        return 0.0
    imax = np.where(stellar.mass == np.max(stellar.mass))[0][0]
    sinkid = stellar.sinkid[imax]-1
    sink = sinks.FindSinks(snap)
    mainsinkmass = sink.mass[sinkid]
    return mainsinkmass
findmainsinkmass = Hamu.Algorithm(_findmainsinkmass)

def plotforsims(labels,plotdiff=False):
    plt.clf()  
    simnames = labels.keys()
    difftxt = ""
    for simname in simnames:
        print("Running for sim", simname)
        sim = hamusims[simname]
        label = labels[simname]
        tcreated = timefuncs.FindTcreatedFirstStar(sim)
        times, sinkmasses = timefuncs.timefunc(sim,findmainsinkmass)
        times -= tcreated
        x = times
        y = sinkmasses
        if plotdiff:
            difftxt = "_diff"
            y = np.diff(sinkmasses) / np.diff(times) / 1000.0 # give in kyr
            x = 0.5*(times[1:] + times[:-1])
        colour = linestyles.Colour(simname)
        line = "-"

        #if "NOFB" in simname:
        #    line = "--"
        plt.gca().plot(x,y,color=colour,label=label,
                linestyle=line,alpha=0.9,                                                                      
                path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()]) 
    plt.legend(fontsize="x-small",frameon=False)
    plt.xlim([0,0.4]) # Myr
    plt.xlabel("Time / Myr")
    plt.ylabel("Main Sink Mass / M$_{\odot}$")
    if plotdiff:
        plt.ylabel("Main Sink Accretion Rate / M$_{\odot} / kyr$")
    plt.savefig("../plots/mainsinkmass"+difftxt+".pdf")




if __name__ == "__main__":
    # Use no feedback runs
    labels = {}
    labels["SEED0_35MSUN_CDMASK_NOFB"] = "Seed0, No Feedback"
    labels["SEED1_35MSUN_CDMASK_NOFB"] = "Seed1, No Feedback"
    labels["SEED2_35MSUN_CDMASK_NOFB"] = "Seed2, No Feedback"
    labels["SEED3_35MSUN_CDMASK_NOFB"] = "Seed3, No Feedback"
    labels["SEED4_35MSUN_CDMASK_NOFB"] = "Seed4, No Feedback"
    #labels["SEED0_35MSUN_CDMASK_WINDUV"] = "Seed0, Wind \& UV"
    #labels["SEED1_35MSUN_CDMASK_WINDUV"] = "Seed1, Wind \& UV"
    #labels["SEED2_35MSUN_CDMASK_WINDUV"] = "Seed2, Wind \& UV"
    #labels["SEED3_35MSUN_CDMASK_WINDUV"] = "Seed3, Wind \& UV"
    #labels["SEED4_35MSUN_CDMASK_WINDUV"] = "Seed4, Wind \& UV"

    plotforsims(labels)
    plotforsims(labels,plotdiff=True)
