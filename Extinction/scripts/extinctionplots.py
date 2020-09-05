"""
Make extinction plots
Sam Geen, May 2020
"""

from startup import *
import rayprof

from pymses.utils import constants as C  

from scipy import interpolate

AVtext = "A$_{\mathrm{V}}$"

def extinctionbyparticle(simname,wholesink):
    # wholesink - if True, use whole sink mass, otherwise list by stellar objects
    sim = hamusims[simname]
    particletracks = {}
    isnaps = {}
    # Hack to get this to run on already processed data
    for snap in sim.Snapshots():
        print("RUNNING FOR SNAP", snap.OutputNumber())
        stellar = stellars.FindStellar(snap.RawData())
        if len(stellar.mass) == 0:
            continue
        sink = sinks.FindSinks(snap)
        # Load columns of each sink
        sinktocolumns = makecolumnsinsinksdict(snap)
        if len(sinktocolumns.keys()) <= 0:
            continue
        # Either list by stellar mass
        if not wholesink:
            particles = stellar.mass
            sinklist = stellar.sinkid
        # Or by sink ID for whole sink
        else:
            particles = sink.id
            sinklist = sink.id
        for particle, sinkid in zip(particles,sinklist):
            if particle not in particletracks:
                particletracks[particle] = []
                isnaps[particle] = []
            particletracks[particle].append(sinktocolumns[sinkid])
            isnaps[particle].append(snap.OutputNumber()-1)
    for particle, track in particletracks.iteritems():
        # Shape will be ntimes x nrays (nrays = 100 by default)
        particletracks[particle] = np.array(track)
    return particletracks, isnaps

def extinctionattimes(simname,times,wholesink):
    sim = hamusims[simname]
    snaptimes = np.array([snap.Time() for snap in sim.Snapshots()])
    snap = sim.Snapshots()[0]
    codetoMyr = snap.RawData().info["unit_time"].express(C.Myr)
    snaptimesMyr = snaptimes*codetoMyr
    if len(times) == 0:
        times = np.linspace(0,snaptimesMyr[-1],1000)
    particletracks, isnaps = extinctionbyparticle(simname,wholesink)
    interptracks = {}
    for particle, track in particletracks.iteritems():
        interptracks[particle] = []
        nrays = track.shape[1]
        for iray in range(0,nrays):
            parttimes = snaptimesMyr[isnaps[particle]]
            interptracks[particle].append(interpolate.interp1d(parttimes,
                                                               track[:,iray]))
    for itime, time in enumerate(times):
        plt.clf()
        mincol = 0.1
        for particle, tracks in interptracks.iteritems():
            columns = []
            if wholesink:
                sink = sinks.FindSinks(snap)
            try:
                for track in tracks:
                    c = track(time)
                    if c < mincol:
                        c = mincol
                    columns.append(c)
                cols = np.sort(columns)
                exts = NHtoAv(cols)
                if not wholesink:
                    mass = particle
                else:
                    mass = sink.mass[sink.id == particle]
                marr = exts*0.0 + mass
                # Make weights that peak at 1
                weights = np.concatenate([np.linspace(0.1,1,len(exts)//2),
                                          np.linspace(1,0.1,len(exts)//2)])
                plt.scatter(marr,exts,c=weights,alpha=0.1,edgecolors='none',cmap="Blues")
            except ValueError:
                pass
        plt.xlabel("Stellar Mass / Msun")
        plt.ylabel("V-band extinction")
        plt.yscale("log")
        plt.xlim([8,120])
        plt.ylim([mincol,1e3])
        plt.text(100,0.2,("%.2f" % time)+" Myr")
        figname = "../plots/extinctions/extinctionattime_"+simname+"_"+str(itime).zfill(5)+".png"
        print("Plotting", figname)
        plt.savefig(figname)

def makecolumnsinsinksdict(snap):
    cols,isinks = rayprof.FindColumnsBySinkInSnap(snap)
    if len(isinks) <= 0:
        return {}
    sinktocolumns = {}
    for columns, isink in zip(cols,isinks):
        sinktocolumns[isink] = columns
    return sinktocolumns
        
def cumulativeextinction(snap,simname,wholesink,dofraction,ylims=None):
    # Plot N(<extinction) vs extinction
    # wholesink : bool - use whole sink or just the stellar objects (massive stars)?
    # If set, do mass visible instead of photon count
    # dofraction : bool - if set, do Nvisible/Ntotal versus just Nvisible
    nrays = rayprof.nprofs # Number of rays cast from each sink
    stellar = stellars.FindStellar(snap.RawData())
    sink = sinks.FindSinks(snap)
    particletocolumns = {}
    if not wholesink:
        if len(stellar.mass) == 0:
            return
        masses = stellar.mass
        sinklist = stellar.sinkid
    # Or by sink ID for whole sink
    else:
        if len(sink.mass) == 0:
            return
        masses = sink.mass
        sinklist = sink.id
    sinktocolumns = makecolumnsinsinksdict(snap)
    # Set up values to plot
    nbins = 200
    extaxis = np.logspace(-1,2,nbins) # It's a pun
    cumuldist = np.zeros((nrays,nbins))
    weights = np.concatenate([np.linspace(0.1,1,nrays//2),
                              np.linspace(1,0.1,nrays//2)])
    cumultot = 0.0
    for mass, sinkid in zip (masses, sinklist):
        try:
            columns = sinktocolumns[sinkid]
        except KeyError:
            print(sinkid, sink.id, sink.mass, sinktocolumns)
            raise KeyError
        cols = np.sort(columns)
        exts = NHtoAv(cols)
        toadd = 1.0
        functext = "N"
        funcunits = ""
        if wholesink:
            toadd = mass
            functext = "M"
            funcunits = " / "+Msolar
        cumultot += toadd
        for icol, ext in enumerate(exts):
           cumuldist[icol,np.where(ext < extaxis)] += toadd
    cumuldist = cumultot - cumuldist
    # Plot
    plt.clf()
    cmap = plt.get_cmap("Blues")
    for icol in range(0,nrays//2):
        plt.plot(extaxis,cumuldist[icol,:],color=cmap(weights[icol]))
    for icol in range(nrays//2,nrays)[::-1]:
        plt.plot(extaxis,cumuldist[icol,:],color=cmap(weights[icol]))
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(ylims)
    #    plt.xlabel(AVtext)
    #    plt.ylabel(functext+"($<$"+AVtext+")"+funcunits)
    time = snap.Time()*snap.RawData().info["unit_time"].express(C.Myr)
    xtxt = extaxis.min()
    ytxt = cumuldist[cumuldist > 0].min()*1.1
    if ylims is not None:
        ytxt = ylims[0]*1.2
    #    plt.text(xtxt,ytxt,("%.2f" % time)+" Myr")
    figname = "../plots/extinctions/cumulextinction"+simname+"_"+str(snap.OutputNumber()).zfill(5)+".png"
    print ("Plotting", figname)
    plt.savefig(figname)
    
def cumulativeextinctionforeachsnap(simname,wholesink,dofraction,ylims=None):
    sim = hamusims[simname]
    for snap in sim.Snapshots():
        cumulativeextinction(snap,simname,wholesink,dofraction,ylims)
           
def extinctionforeachsnap(simname):
    sim = hamusims[simname]
    for snap in sim.Snapshots():
        extinctioninsnap(snap, simname)

def extinctioninsnap(snap, simname):
    stellar = stellars.FindStellar(snap)
    if len(stellar.mass) == 0:
        return
    plt.clf()
    # Load columns of each sink
    sinktocolumns = makecolumnsinsinksdict(snap)
    # Plot stellar objects
    for istellar in range(0,len(stellar.mass)):
        sinkid = stellar.sinkid[istellar]
        mass = stellar.mass[istellar]
        columns = sinktocolumns[sinkid]
        cols = np.sort(columns)
        exts = NHtoAv(cols)
        marr = exts*0.0 + mass
        # Make weights that peak at 1
        weights = np.concatenate([np.linspace(0.1,1,len(exts)//2),
                                  np.linspace(1,0.1,len(exts)//2)])
        plt.scatter(marr,exts,c=weights,alpha=0.1,edgecolors='none')
    plt.xlabel("Stellar Mass / Msun")
    plt.ylabel("V-band extinction")
    plt.yscale("log")
    plt.xlim([8,120])
    plt.ylim([0.1,1e3])
    figname = "../plots/extinctions/extinctioninsnap_"+simname+"_"+str(snap.OutputNumber()).zfill(5)+".png"
    print("Plotting", figname)
    plt.savefig(figname)
    
def extinctionhistogram(simname):
    pass
        
if __name__=="__main__":
    #sim = hamusims["18_LEGO"]
    #timefunc(sim,extinctionhistogram)
    #snap = hamusims["18_LEGO"].Snapshots()[50]
    for simname in ["128_LEGO"]:
        cumulativeextinctionforeachsnap(simname,True,True,[1.0,1e4])
        extinctionattimes(simname,[],False)
        extinctionforeachsnap(simname)
    #extinctionbyparticle("18_LEGO",True)
    #extinctionbyparticle("18_LEGO",False)
