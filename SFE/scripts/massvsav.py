'''
Plot cumulative mass vs Av
Sam Geen, June 2016
'''

from startup import *
import columndensity, freefalltime, linestyles

def MakeForLOS(snap,los,Avcutoff,Avlow=None):
    # Get the column density
    im = columndensity.MakeNHMap(snap.hamusnap,los)
    # Get pixel size
    boxlen = snap.info["boxlen"]
    pixlen = boxlen / im.shape[0] # NOTE: boxlen should be in pc!!
    # Flatten the array and sort it
    dens = np.sort(im.flatten())[::-1]
    # Lower limit
    if Avlow is not None:
        gcm2low = AvtoNH(Avlow)*(mHing/X)
        dens[dens < gcm2low] = 0.0
    # Upper limit
    if Avcutoff is not None:
        gcm2cutoff = AvtoNH(Avcutoff)*(mHing/X)
        dens[dens > gcm2cutoff] = gcm2cutoff
    # Find masses
    umass = (pcincm)**2 / Msuning
    masses = pixlen**2 * dens * umass
    cummass = np.cumsum(masses)
    # Scale im to Av from g/cm^2 via NH
    avs = NHtoAv(dens/(mHing/X))
    # Done!
    return avs, cummass
    
MakeForLOSHamu = Hamu.Algorithm(MakeForLOS)

def PlotForSim(sim,Avcutoff):
    # Load at 2*tff
    tff = freefalltime.Tff(sim)
    # Should be in code units already, so no conversion needed
    snap = sim.FindAtTime(2.0*tff)
    # Make plot lines
    col = linestyles.colour(sim.Name())
    Avlow = 0.1
    av, cmass = MakeForLOSHamu(snap,'x',Avcutoff,Avlow)
    plt.plot(av,cmass,color=col,linestyle = '-')   
    av, cmass = MakeForLOSHamu(snap,'y',Avcutoff,Avlow)
    plt.plot(av,cmass,color=col,linestyle = '--')    
    av, cmass = MakeForLOSHamu(snap,'z',Avcutoff,Avlow)
    plt.plot(av,cmass,color=col,linestyle = ':')

def Run(simnames,Avcutoff=None):
    plt.clf()
    for simname in simnames:
        sim = Hamu.Simulation(simname)
        PlotForSim(sim,Avcutoff)
    # Do axes
    plt.xlabel(r"$A_{v}$ (mag)")
    plt.ylabel(r"Cumulative Mass ("+Msolar+")")
    plt.yscale("log")
    plt.ylim([1,1e6])
    # Do legend
    lines, labels = linestyles.sizelegend()
    plt.legend(lines,labels,fontsize="small",frameon=False,loc="upper right")
    # Save
    cutofftxt = ""
    if Avcutoff is not None:
        cutofftxt = "_"+str(Avcutoff).replace(".","p")
    plt.savefig("../plots/massvsav"+cutofftxt+".pdf")
