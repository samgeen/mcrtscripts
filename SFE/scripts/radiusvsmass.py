'''
Plot radius vs mass for various Av limits
Sam Geen, June 2016
'''

from startup import *
import columndensity, freefalltime, linestyles

def PlotForSim(sim,Akcutoff=None):
    # Load at tff
    tff = freefalltime.Tff(sim)
    # Should be in code units already, so no conversion needed
    snap = sim.FindAtTime(tff)
    # Make plot lines
    col = linestyles.colour(sim.Name())
    def PlotinLOS(snap,los,Akcutoff):
        NHcutoff = None
        if Akcutoff is not None:
            NHcutoff = AktoNH(Akcutoff)
        r, m = columndensity.MassRadiusHamu(snap,los,
                                            NHhigh=NHcutoff)
        plt.plot(r,m,color=col,linestyle = '-')
    PlotinLOS(snap,'x',Akcutoff)
    PlotinLOS(snap,'y',Akcutoff)
    PlotinLOS(snap,'z',Akcutoff)

def Run(simnames,Akcutoff=None):
    plt.clf()
    # Plot for simulations
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        PlotForSim(sim,Akcutoff)
    # Overplot Larson relation M = 100 * r^2
    r = 10.0**np.arange(-1,1,0.2)
    m = 380.0 * r**1.6 # Lombardi 2010, Section 2.4
    plt.plot(r,m,color=r'#666666')
    # Do axes
    plt.xlabel(r"Radius / pc")
    plt.xscale("log")
    plt.ylabel(r"Mass / "+Msolar)
    plt.yscale("log")
    # Do legend
    lines, labels = linestyles.sizelegend()
    plt.legend(lines,labels,fontsize="small",frameon=False,loc="upper left")
    # Save
    cutofftxt = ""
    if Akcutoff is not None:
        cutofftxt = "_"+str(Akcutoff).replace(".","p")
    plt.savefig("../plots/radiusvsmass"+cutofftxt+".pdf")
