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

def plotforsims(labels):
    plt.clf()  
    simnames = labels.keys()
    isim = 0
    for simname in simnames:
        print("Running for sim", simname)
        isim += 1
        #if isim > len(simnames)//2:
        #    linestyles.isim = 0
        sim = hamusims[simname]
        label = labels[simname]
        tcreated = timefuncs.FindTcreatedFirstStar(sim)
        times, sinkmasses = timefuncs.timefunc(sim,findmainsinkmass)
        times -= tcreated
        plt.plot()
        colour = linestyles.Colour(simname)
        line = "-"
        #if "NOFB" in simname:
        #    line = "--"
        plt.gca().plot(times,sinkmasses,color=colour,label=label,
                linestyle=line,alpha=0.9,                                                                      
                path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()]) 
    plt.legend(fontsize="x-small",frameon=False)
    plt.xlim([0,0.4]) # Myr
    plt.xlabel("Time / Myr")
    plt.ylabel("Main Sink Mass / M$_{\odot}$")
    plt.savefig("../plots/mainsinkmass.pdf")




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
